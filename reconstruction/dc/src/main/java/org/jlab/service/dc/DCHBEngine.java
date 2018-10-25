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
import org.jlab.rec.dc.Constants;
import org.jlab.rec.dc.banks.HitReader;
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
import org.jlab.utils.groups.IndexedTable;

// NOTE: Lacks javadoc description
/**
 * @author ziegler
 * @since 08.09.2018 updated by gurjyan
 */
public class DCHBEngine extends DCEngine {

    private AtomicInteger Run = new AtomicInteger(0);
    private double triggerPhase;
    private int newRun = 0;

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

    // NOTE: A detailed javadoc description of the engine itself and what it
    //       does here would greatly benefit anyone joining the project.
    @Override
    public boolean processDataEvent(DataEvent event) {
        // long startTime = 0;
        // setRunConditionsParameters(event);
        if (!event.hasBank("RUN::config")) return true;

        DataBank bank = event.getBank("RUN::config");
        long timeStamp = bank.getLong("timestamp", 0);
        double triggerPhase = 0;

        // Load the constants
        int newRun = bank.getInt("run", 0);
        if (newRun == 0) return true;

        if (Run.get() == 0 || (Run.get() != 0 && Run.get() != newRun)) {
            if (timeStamp == -1) return true;
            // if (debug.get()) startTime = System.currentTimeMillis();
            IndexedTable tabJ = super.getConstantsManager().getConstants(newRun,
                    Constants.TIMEJITTER);
            double period     = tabJ.getDoubleValue("period", 0, 0, 0);
            int phase         = tabJ.getIntValue   ("phase",  0, 0, 0);
            int cycles        = tabJ.getIntValue   ("cycles", 0, 0, 0);

            if (cycles > 0) triggerPhase = period * ((timeStamp + phase) % cycles);

            TableLoader.FillT0Tables(newRun, super.variationName);
            TableLoader.Fill(super.getConstantsManager()
                                  .getConstants(newRun, Constants.TIME2DIST));

            Run.set(newRun);
            if (event.hasBank("MC::Particle") &&
                    this.getEngineConfigString("wireDistort") == null) {
                Constants.setWIREDIST(0);
            }

            // if (debug.get()) {
            //     System.out.println("NEW RUN INIT = " +
            //                        (System.currentTimeMillis() - startTime));
            // }
        }

// === Initial setup ===========================================================
        // Get Field
        Swim dcSwim = new Swim();
        // Init SNR
        Clas12NoiseResult results = new Clas12NoiseResult();
        Clas12NoiseAnalysis noiseAnalysis = new Clas12NoiseAnalysis();
        NoiseReductionParameters parameters =
                new NoiseReductionParameters(2,
                                             Constants.SNR_LEFTSHIFTS,
                                             Constants.SNR_RIGHTSHIFTS);
        ClusterFitter cf           = new ClusterFitter();
        ClusterCleanerUtilities ct = new ClusterCleanerUtilities();
        RecoBankWriter rbc         = new RecoBankWriter();
        HitReader hitRead          = new HitReader();

        hitRead.fetchDCHits(event,
                noiseAnalysis,
                parameters,
                results,
                super.getConstantsManager().getConstants(newRun, Constants.TIME2DIST),
                super.getConstantsManager().getConstants(newRun, Constants.TDCTCUTS),
                super.getConstantsManager().getConstants(newRun, Constants.WIRESTAT),
                dcDetector,
                triggerPhase);

// === Get and setup the hits, clusters, segments and crosses ==================
        // Get the hits
        List<Hit> hits = hitRead.get_DCHits();
        if (hits.isEmpty()) return true;

        // Find the fitted clusters from these hits
        ClusterFinder clusFinder = new ClusterFinder();
        List<FittedCluster> clusters =
                clusFinder.FindHitBasedClusters(hits, ct, cf, dcDetector);
        if (clusters.isEmpty()) return true;

        List<FittedHit> fhits = rbc.createRawHitList(hits);
        rbc.updateHitsListWithClusterInfo(fhits, clusters);

        // Find the segments from the fitted clusters
        SegmentFinder segFinder = new SegmentFinder();
        List<Segment> segments =
                segFinder.get_Segments(clusters, event, dcDetector, false);
        if (segments.isEmpty()) {
            rbc.fillAllHBBanks(event,
                               rbc,
                               fhits,
                               clusters,
                               null,
                               null,
                               null);
            return true;
        }

        // Build a trajectory using six segments
        List<Segment> rmSegs = new ArrayList<>();
        double trkDocOverCellSize;
        for (Segment se : segments) {
            trkDocOverCellSize = 0;
            for (FittedHit fh : se.get_fittedCluster()) {
                trkDocOverCellSize += fh.get_ClusFitDoca() / fh.get_CellSize();
            }
            if (trkDocOverCellSize / se.size() > 1.1) rmSegs.add(se);
        }
        segments.removeAll(rmSegs);

        // Build the crosses using the segments
        CrossMaker crossMake = new CrossMaker();
        List<Cross> crosses = crossMake.find_Crosses(segments, dcDetector);
        if (crosses.isEmpty()) {
            rbc.fillAllHBBanks(event,
                               rbc,
                               fhits,
                               clusters,
                               segments,
                               null,
                               null);
            return true;
        }

        CrossListFinder crossLister = new CrossListFinder();
        CrossList crosslist = crossLister.candCrossLists(crosses,
                false,
                super.getConstantsManager().getConstants(newRun,
                                                         Constants.TIME2DIST),
                dcDetector,
                null,
                dcSwim);

        // Find the list of track candidates using the crosses
        TrackCandListFinder trkcandFinder =
                new TrackCandListFinder(Constants.HITBASE);
        List<Track> trkcands = trkcandFinder.getTrackCands(crosslist,
                                                           dcDetector,
                                                           Swimmer.getTorScale(),
                                                           dcSwim);

        // Process the obtained tracks
        int trkId = 1;
        if (trkcands.size() > 0) {
            // Remove overlaps
            trkcandFinder.removeOverlappingTracks(trkcands);
            for (Track trk : trkcands) {
                // Reset the id
                trk.set_Id(trkId);
                trkcandFinder.matchHits(trk.get_Trajectory(),
                                        trk,
                                        dcDetector,
                                        dcSwim);
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
            if (!c.get_Segment1().isOnTrack)
                crossSegsNotOnTrack.add(c.get_Segment1());
            if (!c.get_Segment2().isOnTrack)
                crossSegsNotOnTrack.add(c.get_Segment2());
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
                    }
                    else {
                        missingSL = r.get(ri).get_Superlayer() - 1;
                    }
                }
            }
            for (int ri = 0; ri < 3; ri++) {
                for (Segment s : crossSegsNotOnTrack) {
                    // if (s.get_Sector() == r.get(ri).get_Sector() &&
                    //         s.get_Region() == r.get(ri).get_Region() &&
                    //         s.associatedCrossId == r.get(ri).associatedCrossId &&
                    //         r.get(ri).associatedCrossId != -1) {
                    //     if (s.get_Superlayer() % 2 == missingSL % 2)
                    //         Segs2Road.add(s);
                    // }
                    if (s.get_Sector() == r.get(ri).get_Sector() &&
                            s.get_Region() == r.get(ri).get_Region() &&
                            s.associatedCrossId == r.get(ri).associatedCrossId &&
                            r.get(ri).associatedCrossId != -1 &&
                            s.get_Superlayer() % 2 == missingSL % 2) {
                        Segs2Road.add(s);
                    }
                }
            }
            if (Segs2Road.size() == 2) {
                Segment pSegment = rf.findRoadMissingSegment(Segs2Road,
                                                             dcDetector,
                                                             r.a);
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
        List<Track> mistrkcands =
                trkcandFinder.getTrackCands(pcrosslist,
                                            dcDetector,
                                            Swimmer.getTorScale(),
                                            dcSwim);

        // Remove overlaps
        if (mistrkcands.size() > 0) {
            trkcandFinder.removeOverlappingTracks(mistrkcands);
            for (Track trk : mistrkcands) {

                // Reset the id
                trk.set_Id(trkId);
                trkcandFinder.matchHits(trk.get_Trajectory(),
                                        trk,
                                        dcDetector,
                                        dcSwim);
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

        // No candidates found. Stop here to save the hits, clusters, segments
        //         and crosses.
        if (trkcands.isEmpty()) {
            rbc.fillAllHBBanks(event,
                               rbc,
                               fhits,
                               clusters,
                               segments,
                               crosses,
                               null);
        }
        else {
            rbc.fillAllHBBanks(event,
                               rbc,
                               fhits,
                               clusters,
                               segments,
                               crosses,
                               trkcands);
        }
        return true;
    }
}
