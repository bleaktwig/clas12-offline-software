package org.jlab.service.dc;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;

import cnuphys.snr.NoiseReductionParameters;
import cnuphys.snr.clas12.Clas12NoiseAnalysis;
import cnuphys.snr.clas12.Clas12NoiseResult;
import org.jlab.clas.swimtools.MagFieldsEngine;
import org.jlab.clas.swimtools.Swim;
import org.jlab.clas.swimtools.Swimmer;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.io.hipo.HipoDataSource;
import org.jlab.io.hipo.HipoDataSync;
import org.jlab.utils.groups.IndexedTable;

import org.jlab.rec.dc.Constants;
import org.jlab.rec.dc.banks.HitReader;
import org.jlab.rec.dc.banks.RecoBankReader;
import org.jlab.rec.dc.banks.RecoBankWriter;
import org.jlab.rec.dc.cluster.ClusterCleanerUtilities;
import org.jlab.rec.dc.cluster.ClusterFinder;
import org.jlab.rec.dc.cluster.ClusterFitter;
import org.jlab.rec.dc.cluster.FittedCluster;
import org.jlab.rec.dc.cross.Cross;
import org.jlab.rec.dc.cross.CrossList;
import org.jlab.rec.dc.cross.CrossListFinder;
import org.jlab.rec.dc.cross.CrossMaker;
import org.jlab.rec.dc.hit.FittedHit;
import org.jlab.rec.dc.hit.Hit;
import org.jlab.rec.dc.segment.Segment;
import org.jlab.rec.dc.segment.SegmentFinder;
import org.jlab.rec.dc.timetodistance.TableLoader;
import org.jlab.rec.dc.track.Track;
import org.jlab.rec.dc.track.TrackCandListFinder;
import org.jlab.rec.dc.trajectory.RoadFinder;
import org.jlab.rec.dc.trajectory.Road;

/**
 * @author ziegler
 * @since 08.09.2018 updated by gurjyan
 */
public class DCHBEngine extends DCEngine {
    private AtomicInteger Run = new AtomicInteger(0);
    private double triggerPhase;
    private int newRun = 0;

    private int eventCounter = 0;

    public DCHBEngine() {
        super("DCHB");
    }

    @Override
    public boolean init() {
        // Load cuts
        Constants.Load();
        super.setStartTimeOption();
        super.LoadTables();
        // newRun = 809;
        // long timeStamp = 371468548086L;
        // if (Run.get() == 0 || (Run.get() != 0 && Run.get() != newRun)) {
        //     IndexedTable tabJ = super.getConstantsManager().getConstants(newRun,
        //             Constants.TIMEJITTER);
        //     double period = tabJ.getDoubleValue("period", 0, 0, 0);
        //     int phase = tabJ.getIntValue("phase", 0, 0, 0);
        //     int cycles = tabJ.getIntValue("cycles", 0, 0, 0);
        //
        //     if (cycles > 0) triggerPhase = period * ((timeStamp + phase) % cycles);
        //
        //     TableLoader.FillT0Tables(newRun, super.variationName);
        //     TableLoader.Fill(super.getConstantsManager().getConstants(newRun, Constants.TIME2DIST));
        //
        //     Run.set(newRun);
        // }
        return true;
    }

    @Override
    public boolean processDataEvent(DataEvent event) {
        int currentEvent = eventCounter;
        eventCounter++;

        // setRunConditionsParameters(event);
        if (!event.hasBank("RUN::config")) {
            return true;
        }

        DataBank bank  = event.getBank("RUN::config");
        long timeStamp = bank.getLong("timestamp", 0);
        double triggerPhase = 0;

        // Load the constants
        int newRun = bank.getInt("run", 0);
        if (newRun == 0) return true;

        if (Run.get() == 0 || (Run.get() != 0 && Run.get() != newRun)) {
            if (timeStamp == -1) return true;
            IndexedTable tabJ = super.getConstantsManager().getConstants(newRun, Constants.TIMEJITTER);
            double period     = tabJ.getDoubleValue("period", 0, 0, 0);
            int phase         = tabJ.getIntValue("phase", 0, 0, 0);
            int cycles        = tabJ.getIntValue("cycles", 0, 0, 0);

            if (cycles > 0) triggerPhase = period * ((timeStamp + phase) % cycles);

            TableLoader.FillT0Tables(newRun, super.variationName);
            TableLoader.Fill(super.getConstantsManager().getConstants(newRun, Constants.TIME2DIST));

            Run.set(newRun);
            if (event.hasBank("MC::Particle") && this.getEngineConfigString("wireDistort") == null) {
                Constants.setWIREDIST(0);
            }
        }

        // get Field
        Swim dcSwim = new Swim();

        // init SNR
        Clas12NoiseResult results = new Clas12NoiseResult();
        Clas12NoiseAnalysis noiseAnalysis = new Clas12NoiseAnalysis();
        NoiseReductionParameters parameters =
                new NoiseReductionParameters(2, Constants.SNR_LEFTSHIFTS, Constants.SNR_RIGHTSHIFTS);

        ClusterFitter cf = new ClusterFitter();
        ClusterCleanerUtilities ct = new ClusterCleanerUtilities();
        RecoBankWriter rbc = new RecoBankWriter();
        HitReader hitRead = new HitReader();
        hitRead.fetch_DCHits(event,
                noiseAnalysis,
                parameters,
                results,
                super.getConstantsManager().getConstants(newRun, Constants.TIME2DIST),
                super.getConstantsManager().getConstants(newRun, Constants.TDCTCUTS),
                super.getConstantsManager().getConstants(newRun, Constants.WIRESTAT),
                dcDetector,
                triggerPhase);

        // Hits
        List<Hit> hits = hitRead.get_DCHits();
        if (hits.isEmpty()) return true;

        // Clusters
        ClusterFinder clusFinder = new ClusterFinder();
        List<FittedCluster> clusters = clusFinder.FindHitBasedClusters(hits, ct, cf, dcDetector);
        if (clusters.isEmpty()) return true;
        List<FittedHit> fhits = rbc.createRawHitList(hits);
        rbc.updateHitsListWithClusterInfo(fhits, clusters);

        // Segments
        SegmentFinder segFinder = new SegmentFinder();
        List<Segment> segments = segFinder.get_Segments(clusters, event, dcDetector, false);

        // Build trajectory using segments
        if (segments.isEmpty()) {
            rbc.fillAllBanks(event, rbc, fhits, clusters, null, null, null, false);
            return true;
        }
        List<Segment> rmSegs = new ArrayList<>();

        // Clean up hit-based segments
        double trkDocOverCellSize;
        for (Segment se : segments) {
            trkDocOverCellSize = 0;
            for (FittedHit fh : se.get_fittedCluster()) {
                trkDocOverCellSize += fh.get_ClusFitDoca() / fh.get_CellSize();
            }
            if (trkDocOverCellSize / se.size() > 1.1) {
                rmSegs.add(se);
            }
        }
        segments.removeAll(rmSegs);

        // Crosses
        CrossMaker crossMake = new CrossMaker();
        List<Cross> crosses = crossMake.find_Crosses(segments, dcDetector);

        // if (currentEvent == 0) {
        //     System.out.println("\n\n\n\nDCHB DATA:\n");
        //     RecoBankReader.printSample(crosses.get(0));
        // }

        // // TODO: Solve the chi2 inconsistency
        // if (currentEvent >= 3 && currentEvent <= 5) {
        //     System.out.println("[DCHB." + currentEvent + "] - 2 - cluster "
        //                     + crosses.get(0).get_Segment1().get_fittedCluster().get_Id()
        //                     + "'s fit chi2: "
        //                     + crosses.get(0).get_Segment1().get_fittedCluster().get_fitProb());
        // }

        if (crosses.isEmpty()) {
            rbc.fillAllBanks(event, rbc, fhits, clusters, segments, null, null, false);
            return true;
        }
        else {
            rbc.fillAllBanks(event, rbc, fhits, clusters, segments, crosses, null, false);
        }

        /* 17 */
        CrossListFinder crossLister = new CrossListFinder();

        CrossList crosslist = crossLister.candCrossLists(crosses,
                false,
                super.getConstantsManager().getConstants(newRun, Constants.TIME2DIST),
                dcDetector,
                null,
                dcSwim);

        if (currentEvent == -1) {
            /* 18 */
            //6) find the list of  track candidates
            TrackCandListFinder trkcandFinder = new TrackCandListFinder(Constants.HITBASE);
            List<Track> trkcands = trkcandFinder.getTrackCands(crosslist,
                    dcDetector,
                    Swimmer.getTorScale(),
                    dcSwim);
            /* 19 */

            // track found
            int trkId = 1;
            if (trkcands.size() > 0) {
                // remove overlaps
                trkcandFinder.removeOverlappingTracks(trkcands);
                for (Track trk : trkcands) {
                    // reset the id
                    trk.set_Id(trkId);
                    trkcandFinder.matchHits(trk.get_Trajectory(), trk, dcDetector, dcSwim);
                    for (Cross c : trk) {
                        c.get_Segment1().isOnTrack = true;
                        c.get_Segment2().isOnTrack = true;

                        for (FittedHit h1 : c.get_Segment1()) {
                            h1.set_AssociatedHBTrackID(trk.get_Id());
                        }
                        for (FittedHit h2 : c.get_Segment2()) {
                            h2.set_AssociatedHBTrackID(trk.get_Id());
                        }
                    }
                    trkId++;
                }
            }
            List<Segment> crossSegsNotOnTrack = new ArrayList<>();
            List<Segment> psegments = new ArrayList<>();

            for (Cross c : crosses) {
                if (!c.get_Segment1().isOnTrack) crossSegsNotOnTrack.add(c.get_Segment1());
                if (!c.get_Segment2().isOnTrack) crossSegsNotOnTrack.add(c.get_Segment2());
            }

            RoadFinder rf = new RoadFinder();
            List<Road> allRoads = rf.findRoads(segments, dcDetector);
            List<Segment> Segs2Road = new ArrayList<>();
            for (Road r : allRoads) {
                Segs2Road.clear();
                int missingSL = -1;
                for (int ri = 0; ri < 3; ri++) {
                    if (r.get(ri).associatedCrossId == -1) {
                        if (r.get(ri).get_Superlayer() % 2 == 1) {
                            missingSL = r.get(ri).get_Superlayer() + 1;
                        } else {
                            missingSL = r.get(ri).get_Superlayer() - 1;
                        }
                    }
                }
                for (int ri = 0; ri < 3; ri++) {
                    for (Segment s : crossSegsNotOnTrack) {
                        if (s.get_Sector() == r.get(ri).get_Sector()
                                && s.get_Region() == r.get(ri).get_Region()
                                && s.associatedCrossId == r.get(ri).associatedCrossId
                                && r.get(ri).associatedCrossId != -1) {

                            if (s.get_Superlayer() % 2 == missingSL % 2) Segs2Road.add(s);
                        }
                    }
                }
                if (Segs2Road.size() == 2) {
                    Segment pSegment = rf.findRoadMissingSegment(Segs2Road, dcDetector, r.a);
                    if (pSegment != null) psegments.add(pSegment);
                }
            }

            segments.addAll(psegments);
            List<Cross> pcrosses = crossMake.find_Crosses(segments, dcDetector);
            CrossList pcrosslist = crossLister.candCrossLists(pcrosses,
                    false,
                    super.getConstantsManager().getConstants(newRun, Constants.TIME2DIST),
                    dcDetector,
                    null,
                    dcSwim);
            List<Track> mistrkcands = trkcandFinder.getTrackCands(pcrosslist,
                    dcDetector,
                    Swimmer.getTorScale(),
                    dcSwim);

            // remove overlaps
            if (mistrkcands.size() > 0) {
                trkcandFinder.removeOverlappingTracks(mistrkcands);
                for (Track trk : mistrkcands) {

                    // reset the id
                    trk.set_Id(trkId);
                    trkcandFinder.matchHits(trk.get_Trajectory(), trk, dcDetector, dcSwim);
                    for (Cross c : trk) {
                        for (FittedHit h1 : c.get_Segment1()) {
                            h1.set_AssociatedHBTrackID(trk.get_Id());
                        }
                        for (FittedHit h2 : c.get_Segment2()) {
                            h2.set_AssociatedHBTrackID(trk.get_Id());
                        }
                    }
                    trkId++;
                }
            }
            trkcands.addAll(mistrkcands);

            // no candidate found, stop here and save the hits,
            // the clusters, the segments, the crosses
            if (trkcands.isEmpty()) {
                rbc.fillAllBanks(event, rbc, fhits, clusters, segments, crosses, null, false);
                return true;
            }
            rbc.fillAllBanks(event, rbc, fhits, clusters, segments, crosses, trkcands, false);
        }
        return true;
    }
}
