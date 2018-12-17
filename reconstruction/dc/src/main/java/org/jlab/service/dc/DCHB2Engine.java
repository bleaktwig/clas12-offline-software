package org.jlab.service.dc;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;

// import cnuphys.snr.NoiseReductionParameters;
// import cnuphys.snr.clas12.Clas12NoiseAnalysis;
// import cnuphys.snr.clas12.Clas12NoiseResult;
// import org.jlab.clas.swimtools.MagFieldsEngine;
import org.jlab.clas.swimtools.Swim;
import org.jlab.clas.swimtools.Swimmer;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
// import org.jlab.io.hipo.HipoDataSource;
// import org.jlab.io.hipo.HipoDataSync;
// import org.jlab.utils.groups.IndexedTable;

import org.jlab.rec.dc.Constants;
// import org.jlab.rec.dc.banks.HitReader;
import org.jlab.rec.dc.banks.RecoBankReader;
import org.jlab.rec.dc.banks.RecoBankWriter;
// import org.jlab.rec.dc.cluster.ClusterCleanerUtilities;
// import org.jlab.rec.dc.cluster.ClusterFinder;
// import org.jlab.rec.dc.cluster.ClusterFitter;
import org.jlab.rec.dc.cluster.FittedCluster;
import org.jlab.rec.dc.cross.Cross;
import org.jlab.rec.dc.cross.CrossList;
import org.jlab.rec.dc.cross.CrossListFinder;
import org.jlab.rec.dc.cross.CrossMaker;
import org.jlab.rec.dc.hit.FittedHit;
// import org.jlab.rec.dc.hit.Hit;
import org.jlab.rec.dc.segment.Segment;
// import org.jlab.rec.dc.segment.SegmentFinder;
import org.jlab.rec.dc.timetodistance.TableLoader;
import org.jlab.rec.dc.track.Track;
import org.jlab.rec.dc.track.TrackCandListFinder;
import org.jlab.rec.dc.trajectory.RoadFinder;
import org.jlab.rec.dc.trajectory.Road;

/**
 * @author ziegler
 * @since 08.09.2018 updated by gurjyan
 */
public class DCHB2Engine extends DCEngine {
    private AtomicInteger Run = new AtomicInteger(0);
    private double triggerPhase;
    private int newRun = 0;

    private int eventCounter = 0;

    public DCHB2Engine() {
        super("DCHB2");
    }

    @Override
    public boolean init() {
        // Load cuts
        Constants.Load();
        super.setStartTimeOption();
        super.LoadTables();
        return true;
    }

    @Override
    public boolean processDataEvent(DataEvent event) {
        int currentEvent = eventCounter;
        eventCounter++;
        if (currentEvent != 21) return true;

        // === INITIAL CHECKUP =========================================================
        if (!event.hasBank("RUN::config")) return true;
        DataBank headerBank = event.getBank("RUN::config");

        int newRun = headerBank.getInt("run", 0);
        if (newRun == 0) return true;

        Swim dcSwim = new Swim();
        RecoBankWriter rbw = new RecoBankWriter();
        RecoBankReader rbr = new RecoBankReader();

        CrossMaker crossMake = new CrossMaker();
        CrossListFinder crossLister = new CrossListFinder();
        TrackCandListFinder trkCandFinder = new TrackCandListFinder(Constants.HITBASE);

        // === GET OBJECTS =============================================================
        if (!event.hasBank("HitBasedTrkg::HBCrosses")) return true;

        DataBank hiBank = event.getBank("HitBasedTrkg::HBHits");
        DataBank clBank = event.getBank("HitBasedTrkg::HBClusters");
        DataBank seBank = event.getBank("HitBasedTrkg::HBSegments");
        DataBank crBank = event.getBank("HitBasedTrkg::HBCrosses");

        // Pull the hits, clusters, segments and crosses from their banks
        if (hiBank.rows() == 0) return true;
        if (clBank.rows() == 0) return true;
        if (seBank.rows() == 0) return true;
        if (crBank.rows() == 0) return true;

        List<FittedHit> hits         = new ArrayList();
        List<FittedCluster> clusters = new ArrayList();
        List<Segment> segments       = new ArrayList();
        List<Cross> crosses          = new ArrayList();

        for (int h = 0; h < hiBank.rows(); ++h) hits.add(rbr.getHit(hiBank, h));
        for (int c = 0; c < clBank.rows(); ++c) clusters.add(rbr.getCluster(clBank, hits, c));
        for (int s = 0; s < seBank.rows(); ++s) segments.add(rbr.getSegment(seBank, clusters, s));
        for (int c = 0; c < crBank.rows(); ++c) crosses.add(rbr.getCross(crBank, segments, c));

        // Pull the tracks from the banks if available
        List<Track> trkCands = new ArrayList();
        DataBank trBank = event.getBank("HitBasedTrkg::HBTracks");
        if (trBank.rows() > 0) {
            for (int t = 0; t < trBank.rows(); ++t) trkCands.add(rbr.getHBTrack(trBank, crosses, t));
        }
        else System.out.println("YES!");

        // === RUN DCHB2 ===============================================================
        // System.out.println("[DCHB2] 00 track candidates size: " + trkCands.size());
        int trkId = 1;
        if (trkCands.size() > 0) {
            for (Track trk : trkCands) {
                trk.set_Id(trkId); // Reset the id
                // System.out.println("[DCHB2] 01 matching hits");
                trkCandFinder.matchHits(trk.get_Trajectory(), trk, dcDetector, dcSwim);
                for (Cross c : trk) {
                    c.get_Segment1().isOnTrack = true;
                    c.get_Segment2().isOnTrack = true;

                    for (FittedHit h1 : c.get_Segment1()) h1.set_AssociatedHBTrackID(trk.get_Id());
                    for (FittedHit h2 : c.get_Segment2()) h2.set_AssociatedHBTrackID(trk.get_Id());
                }
                trkId++;
            }
        }

        List<Segment> crossSegsNotOnTrack = new ArrayList<>();
        List<Segment> psegments           = new ArrayList<>();

        for (Cross c : crosses) {
            if (!c.get_Segment1().isOnTrack) {
                crossSegsNotOnTrack.add(c.get_Segment1());
                // System.out.println("[DCHB2] 01A Segment1 not on track!");
            }
            if (!c.get_Segment2().isOnTrack) {
                crossSegsNotOnTrack.add(c.get_Segment2());
                // System.out.println("[DCHB2] 01B Segment2 not on track!");
            }
        }

        RoadFinder rf = new RoadFinder();
        // System.out.println("[DCHB2] 02 Finding roads");
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

        // System.out.println("[DCHB2] 03 Roads found, making segments");
        segments.addAll(psegments);
        // System.out.println("[DCHB2] 03A psegments size: " + psegments.size());
        // System.out.println("[DCHB2] 03B segments size:  " + segments.size());
        // RecoBankReader.printSampleSegment(segments.get(0));
        // System.out.println("[DCHB2] 04 Making crosses");
        List<Cross> pcrosses = crossMake.find_Crosses(segments, dcDetector);

        // System.out.println("[DCHB2] 04A prcrosses size: " + pcrosses.size());
        // RecoBankReader.printSample(pcrosses.get(0));
        // System.out.println("[DCHB2] 05 Making crosslists");
        CrossList pcrosslist = crossLister.candCrossLists(pcrosses,
                false,
                super.getConstantsManager().getConstants(newRun, Constants.TIME2DIST),
                dcDetector,
                null,
                dcSwim);
        // System.out.println("[DCHB2] 06 Finding mistrkcands");
        List<Track> mistrkcands = trkCandFinder.getTrackCands(pcrosslist,
                dcDetector,
                Swimmer.getTorScale(),
                dcSwim);

        // System.out.println("[DCHB2] 07 # of mistrkcands before removing overlapping tracks: "
                           // + mistrkcands.size());
        // remove overlaps
        if (mistrkcands.size() > 0) {
            trkCandFinder.removeOverlappingTracks(mistrkcands);
            for (Track trk : mistrkcands) {

                // reset the id
                trk.set_Id(trkId);
                trkCandFinder.matchHits(trk.get_Trajectory(), trk, dcDetector, dcSwim);
                for (Cross c : trk) {
                    for (FittedHit h1 : c.get_Segment1()) h1.set_AssociatedHBTrackID(trk.get_Id());
                    for (FittedHit h2 : c.get_Segment2()) h2.set_AssociatedHBTrackID(trk.get_Id());
                }
                trkId++;
            }
        }

        System.out.println("[DCHB2] # of mistrkcands after removing overlapping tracks: "
                           + mistrkcands.size());
        trkCands.addAll(mistrkcands);

        // TODO: Some changes were made to at least the hits, segments and crosses. Remove these
        //       from the banks and write them again.

        System.out.println("[DCHB2] final # of tracks: " + trkCands.size() + "\n");
        if (!trkCands.isEmpty()) rbw.fillHBTracksBanks(event, rbw, trkCands);

        return true;
    }
}
