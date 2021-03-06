 package org.jlab.service.dc;

import cnuphys.snr.NoiseReductionParameters;
import cnuphys.snr.clas12.Clas12NoiseAnalysis;
import cnuphys.snr.clas12.Clas12NoiseResult;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;

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
import org.jlab.rec.dc.track.BFieldInterpolator;
import org.jlab.rec.dc.track.Track;
import org.jlab.rec.dc.track.TrackCandListFinder;
import org.jlab.rec.dc.trajectory.RoadFinder;
import org.jlab.rec.dc.trajectory.Road;
import org.jlab.utils.groups.IndexedTable;

/**
 * @author zigler
 * @since 08.09.2018 updated by gurjyan
 */
public class DCHBEngine extends DCEngine {

    private AtomicInteger Run = new AtomicInteger(0);
    private double triggerPhase;
    private int newRun = 0;

    // Define if cluster and track finding is done sequentially or in parallel:
    //     false: sequential
    //     true:  parallel
    private boolean clusterFindingMode = false;
    private boolean trackFindingMode   = false;

    private BFieldInterpolator[] bField;
    private int eventCtr = 0;
    private boolean firstEvent = true;

    public DCHBEngine() {
        super("DCHB");
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
        if (!event.hasBank("RUN::config")) return true;

        // eventCtr++;
        // boolean runEvent = false;
        // NOTE: Singular covariance matrix events:
        // {
        //     if (eventCtr ==    1) runEvent = true;
        //     if (eventCtr ==    2) runEvent = true;
        //     if (eventCtr ==    3) runEvent = true;
        //
        //     if (eventCtr ==   50) runEvent = true;
        //     if (eventCtr ==  356) runEvent = true;
        //     if (eventCtr ==  695) runEvent = true;
        //     if (eventCtr ==  923) runEvent = true;
        //     if (eventCtr == 1458) runEvent = true;
        //     if (eventCtr == 2679) runEvent = true;
        //     if (eventCtr == 3242) runEvent = true;
        //     if (eventCtr == 4166) runEvent = true;
        //     if (eventCtr == 5221) runEvent = true;
        //     if (eventCtr == 6638) runEvent = true;
        //     if (eventCtr == 7320) runEvent = true;
        //     if (eventCtr == 7689) runEvent = true;
        // }
        // NOTE: One track events:
        // {
        //     if (eventCtr ==    1) runEvent = true;
        //     if (eventCtr ==    2) runEvent = true;
        //     if (eventCtr ==    3) runEvent = true;

        //     if (eventCtr ==    5) runEvent = true;
        //     if (eventCtr ==    6) runEvent = true;
        //     if (eventCtr ==   12) runEvent = true;
        //     if (eventCtr ==   17) runEvent = true;
        //     if (eventCtr ==   25) runEvent = true;
        //     if (eventCtr ==   30) runEvent = true;
        //
        //     if (eventCtr ==   35) runEvent = true;
        //     if (eventCtr ==   36) runEvent = true;
        //     if (eventCtr ==   38) runEvent = true;
        //     if (eventCtr ==   40) runEvent = true;
        //     if (eventCtr ==   48) runEvent = true;
        //     if (eventCtr ==   52) runEvent = true;
        //
        //     if (eventCtr ==   55) runEvent = true;
        //     if (eventCtr ==   57) runEvent = true;
        //     if (eventCtr ==   59) runEvent = true;
        //     if (eventCtr ==   60) runEvent = true;
        //     if (eventCtr ==   63) runEvent = true;
        //     if (eventCtr ==   67) runEvent = true;
        //
        //     if (eventCtr ==   69) runEvent = true;
        //     if (eventCtr ==   71) runEvent = true;
        //     if (eventCtr ==   80) runEvent = true;
        //     if (eventCtr ==   82) runEvent = true;
        //     if (eventCtr ==   83) runEvent = true;
        //     if (eventCtr ==   88) runEvent = true;
        //
        //     if (eventCtr ==   89) runEvent = true;
        //     if (eventCtr ==   92) runEvent = true;
        //     if (eventCtr ==   95) runEvent = true;
        //     if (eventCtr ==  103) runEvent = true;
        //     if (eventCtr ==  113) runEvent = true;
        //     if (eventCtr ==  119) runEvent = true;
        //
        //     if (eventCtr ==  123) runEvent = true;
        //     if (eventCtr ==  124) runEvent = true;
        //     if (eventCtr ==  130) runEvent = true;
        //     if (eventCtr ==  134) runEvent = true;
        //     if (eventCtr ==  138) runEvent = true;
        //     if (eventCtr ==  139) runEvent = true;
        // }
        // if (!runEvent) return true;
        // System.out.printf("# event #%5d\n", eventCtr);

        // long start = System.currentTimeMillis();

        DataBank bank = event.getBank("RUN::config");
        long timeStamp = bank.getLong("timestamp", 0);
        double triggerPhase = 0;

        // Load the constants
        int newRun = bank.getInt("run", 0);
        if (newRun == 0) return true;

        if (Run.get() == 0 || (Run.get() != 0 && Run.get() != newRun)) {
            if (timeStamp == -1) return true;

            IndexedTable tabJ = super.getConstantsManager().getConstants(newRun, Constants.TIMEJITTER);
            double period = tabJ.getDoubleValue("period", 0, 0, 0);
            int phase     = tabJ.getIntValue   ("phase",  0, 0, 0);
            int cycles    = tabJ.getIntValue   ("cycles", 0, 0, 0);

            if (cycles > 0) triggerPhase = period * ((timeStamp + phase) % cycles);

            TableLoader.FillT0Tables(newRun, super.variationName);
            TableLoader.Fill(super.getConstantsManager().getConstants(newRun, Constants.TIME2DIST));

            Run.set(newRun);
            if (event.hasBank("MC::Particle") && this.getEngineConfigString("wireDistort") == null) {
                Constants.setWIREDIST(0);
            }
        }

        Swim dcSwim = new Swim();

        if (firstEvent) {
            // Create bfield interpolators
            double[] min = new double[] {-200., -200., 200.};
            double[] max = new double[] { 200.,  200., 600.};
            double[] ss  = new double[] {  10.,   10.,  10.};

            bField = new BFieldInterpolator[6];
            for (int i = 1; i <= 6; ++i) bField[i-1] = new BFieldInterpolator(dcSwim, i, min, max, ss);

            firstEvent = false;
        }

        // if (eventCtr == 5) {
        //     float[]  fb = new float[3];
        //     double[] rb = new double[3];
        //
        //     for (double x = -197.5; x <= 197.5; x +=   5.) {
        //         for (double y = -197.5; y <= 197.5; y +=   5.) {
        //             for (double z = 202.5; z <= 597.5; z +=   5.) {
        //                 dcSwim.Bfield(2, x, y, z, fb);
        //                 // for (int i = 0; i < fb.length; ++i) rb[i] = (double) fb[i];
        //                 rb = bField[1].getB(x, y, z);
        //                 System.out.printf("%20.12e %20.12e %20.12e\n", rb[0], rb[1], rb[2]);
        //             }
        //         }
        //     }
        // }

        // if (!firstEvent) return true;

        Clas12NoiseResult results = new Clas12NoiseResult();
        Clas12NoiseAnalysis noiseAnalysis = new Clas12NoiseAnalysis();
        NoiseReductionParameters parameters =
                new NoiseReductionParameters(2, Constants.SNR_LEFTSHIFTS, Constants.SNR_RIGHTSHIFTS);

        RecoBankWriter rbc = new RecoBankWriter();

        // Read hits
        HitReader hitRead = new HitReader();

        hitRead.fetchDCHits(event, noiseAnalysis, parameters, results,
                super.getConstantsManager().getConstants(newRun, Constants.TDCTCUTS),
                super.getConstantsManager().getConstants(newRun, Constants.WIRESTAT),
                dcDetector, triggerPhase);

        List<Hit> hits = hitRead.get_DCHits();
        if (hits.isEmpty()) return true;

        // Find the clusters from these hits
        ClusterFitter cf = new ClusterFitter();
        ClusterCleanerUtilities ct = new ClusterCleanerUtilities();

        ClusterFinder clusFinder = new ClusterFinder();
        List<FittedCluster> clusters = clusFinder.FindHitBasedClusters(hits, ct, cf, dcDetector,
                clusterFindingMode);
        if (clusters.isEmpty()) return true;

        // Update hits with cluster information
        List<FittedHit> fhits = rbc.createRawHitList(hits);
        rbc.updateHitsListWithClusterInfo(fhits, clusters);

        // Find the segments from the fitted clusters
        SegmentFinder segFinder = new SegmentFinder();
        List<Segment> segments = segFinder.get_Segments(clusters, event, dcDetector, false);
        if (segments.isEmpty()) {
            rbc.fillAllHBBanks(event, rbc, fhits, clusters, null, null, null);
            return true;
        }
        List<Segment> rmSegs = new ArrayList<>();

        double trkDocOverCellSize;
        for (Segment se : segments) {
            trkDocOverCellSize = 0;
            for (FittedHit fh : se.get_fittedCluster())
                trkDocOverCellSize += fh.get_ClusFitDoca() / fh.get_CellSize();
            if (trkDocOverCellSize / se.size() > 1.1) rmSegs.add(se);
        }
        segments.removeAll(rmSegs);

        // Make the crosses from the segments
        CrossMaker crossMake = new CrossMaker();
        List<Cross> crosses = crossMake.find_Crosses(segments, dcDetector);
        if (crosses.isEmpty()) {
            rbc.fillAllHBBanks(event, rbc, fhits, clusters, segments, null, null);
            return true;
        }

        CrossListFinder crossLister = new CrossListFinder();
        CrossList crosslist = crossLister.candCrossLists(crosses, false,
                super.getConstantsManager().getConstants(newRun, Constants.TIME2DIST),
                dcDetector, null, dcSwim);

        // if (eventCtr <= 4) return true;

        // Find the list of track candidates
        TrackCandListFinder trkcandFinder = new TrackCandListFinder(Constants.HITBASE);
        List<Track> trkcands = trkcandFinder.getTrackCands(crosslist, dcDetector,
                Swimmer.getTorScale(), dcSwim, bField, trackFindingMode);

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

                    for (FittedHit h1 : c.get_Segment1()) h1.set_AssociatedHBTrackID(trk.get_Id());
                    for (FittedHit h2 : c.get_Segment2()) h2.set_AssociatedHBTrackID(trk.get_Id());
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
                            && r.get(ri).associatedCrossId != -1
                            && s.get_Superlayer() % 2 == missingSL % 2)
                        Segs2Road.add(s);
                }
            }
            if (Segs2Road.size() == 2) {
                Segment pSegment = rf.findRoadMissingSegment(Segs2Road, dcDetector, r.a);
                if (pSegment != null) psegments.add(pSegment);
            }
        }
        segments.addAll(psegments);
        List<Cross> pcrosses = crossMake.find_Crosses(segments, dcDetector);
        CrossList pcrosslist = crossLister.candCrossLists(pcrosses, false,
                super.getConstantsManager().getConstants(newRun, Constants.TIME2DIST),
                dcDetector, null, dcSwim);
        List<Track> mistrkcands = trkcandFinder.getTrackCands(pcrosslist, dcDetector,
                Swimmer.getTorScale(), dcSwim, bField, trackFindingMode);

        // remove overlaps
        if (mistrkcands.size() > 0) {
            trkcandFinder.removeOverlappingTracks(mistrkcands);
            for (Track trk : mistrkcands) {

                // reset the id
                trk.set_Id(trkId);
                trkcandFinder.matchHits(trk.get_Trajectory(), trk, dcDetector, dcSwim);
                for (Cross c : trk) {
                    for (FittedHit h1 : c.get_Segment1()) h1.set_AssociatedHBTrackID(trk.get_Id());
                    for (FittedHit h2 : c.get_Segment2()) h2.set_AssociatedHBTrackID(trk.get_Id());
                }
                trkId++;
            }
        }
        trkcands.addAll(mistrkcands);

        // No candidate found, stop here and save the hits, the clusters, the segments, the crosses
        if (trkcands.isEmpty())
            rbc.fillAllHBBanks(event, rbc, fhits, clusters, segments, crosses, null);
        else
            rbc.fillAllHBBanks(event, rbc, fhits, clusters, segments, crosses, trkcands);

        // long end = System.currentTimeMillis();
        // System.out.printf("%d\n", end-start);
        // for (Track trk : trkcands) {
        //     System.out.printf("DCHB TRACK FOUND!!!\n");
        // }

        return true;
    }

    public static void main(String[] args) {

        String inputFile = "/Users/ziegler/Desktop/Work/validation/infiles/straight.hipo";
        MagFieldsEngine enf = new MagFieldsEngine();
        enf.init();

        DCHBEngine en = new DCHBEngine();
        en.init();

        DCTBEngine en2 = new DCTBEngine();
        en2.init();

        int counter = 0;

        HipoDataSource reader = new HipoDataSource();
        reader.open(inputFile);

        HipoDataSync writer = new HipoDataSync();

        String outputFile = "/Users/ziegler/Desktop/Work/Files/test.hipo";

        writer.open(outputFile);
        long t1 = 0;
        while (reader.hasEvent()) {

            counter++;
            System.out.println("************************************************************* ");
            DataEvent event = reader.getNextEvent();
            if (counter > 0) {
                t1 = System.currentTimeMillis();
            }
            enf.processDataEvent(event);
            en.processDataEvent(event);

            // Processing TB
            en2.processDataEvent(event);
            writer.writeEvent(event);
            System.out.println("PROCESSED  EVENT " + event.getBank("RUN::config").getInt("event", 0));

            if(counter>40)
                break;
        }
        writer.close();
        double t = System.currentTimeMillis() - t1;
        System.out.println(t1 + " TOTAL  PROCESSING TIME = " + (t / (float) counter));
    }

}
