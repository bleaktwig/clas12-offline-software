package org.jlab.service.dc;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;

import cnuphys.snr.NoiseReductionParameters;
import cnuphys.snr.clas12.Clas12NoiseAnalysis;
import cnuphys.snr.clas12.Clas12NoiseResult;
import org.jlab.clas.swimtools.Swim;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.utils.groups.IndexedTable;

import org.jlab.rec.dc.Constants;
import org.jlab.rec.dc.banks.HitReader;
import org.jlab.rec.dc.banks.RecoBankWriter;
import org.jlab.rec.dc.timetodistance.TableLoader;
import org.jlab.rec.dc.hit.Hit;
import org.jlab.rec.dc.hit.FittedHit;
import org.jlab.rec.dc.cluster.FittedCluster;
import org.jlab.rec.dc.cluster.ClusterCleanerUtilities;
import org.jlab.rec.dc.cluster.ClusterFinder;
import org.jlab.rec.dc.cluster.ClusterFitter;
import org.jlab.rec.dc.segment.Segment;
import org.jlab.rec.dc.segment.SegmentFinder;
import org.jlab.rec.dc.cross.Cross;
import org.jlab.rec.dc.cross.CrossMaker;

/**
 * @author ziegler
 * @since 08.09.2018 updated by gurjyan
 */
public class DCHB1Engine extends DCEngine {
    private AtomicInteger Run = new AtomicInteger(0);
    private double triggerPhase;
    private int newRun = 0;

    private int eventCounter = 0;

    public DCHB1Engine() {
        super("DCHB1");
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
        if (currentEvent != 42) return true;
        System.out.println("DCHB1 RUNNING on run " + currentEvent);

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

        // Get Field
        Swim dcSwim = new Swim();

        // Init SNR
        Clas12NoiseResult results = new Clas12NoiseResult();
        Clas12NoiseAnalysis noiseAnalysis = new Clas12NoiseAnalysis();
        NoiseReductionParameters parameters =
                new NoiseReductionParameters(2, Constants.SNR_LEFTSHIFTS, Constants.SNR_RIGHTSHIFTS);

        ClusterFitter cf = new ClusterFitter();
        ClusterCleanerUtilities ct = new ClusterCleanerUtilities();
        RecoBankWriter rbw = new RecoBankWriter();

        // Read the hits
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

        List<Hit> hits = hitRead.get_DCHits();
        if (hits.isEmpty()) return true;

        // Find the clusters
        ClusterFinder clusFinder = new ClusterFinder();
        List<FittedCluster> clusters = clusFinder.FindHitBasedClusters(hits, ct, cf, dcDetector);
        if (clusters.isEmpty()) return true;
        List<FittedHit> fhits = rbw.createRawHitList(hits);
        rbw.updateHitsListWithClusterInfo(fhits, clusters);

        // Form the segments from the clusters
        SegmentFinder segFinder = new SegmentFinder();
        List<Segment> segments = segFinder.get_Segments(clusters, event, dcDetector, false);

        // Build trajectory using segments
        if (segments.isEmpty()) {
            rbw.fillAllHBBanks(event, rbw, fhits, clusters, null, null, null);
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

        if (crosses.isEmpty()) {
            rbw.fillAllHBBanks(event, rbw, fhits, clusters, segments, null, null);
        }
        else {
            rbw.fillAllHBBanks(event, rbw, fhits, clusters, segments, crosses, null);
        }

        return true;
    }
}
