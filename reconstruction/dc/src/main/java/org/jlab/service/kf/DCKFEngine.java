package org.jlab.service.kf;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Optional;

import org.jlab.geom.prim.Point3D;
import org.jlab.geom.base.ConstantProvider;
import org.jlab.detector.base.DetectorType;
import org.jlab.detector.base.GeometryFactory;
import org.jlab.detector.geant4.v2.DCGeant4Factory;
import org.jlab.clas.swimtools.Swim;
import org.jlab.clas.swimtools.Swimmer;
import org.jlab.clas.reco.ReconstructionEngine;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;

import org.jlab.rec.dc.Constants;
import org.jlab.rec.dc.banks.RecoBankReader;
import org.jlab.rec.dc.banks.RecoBankWriter;
import org.jlab.rec.dc.hit.FittedHit;
import org.jlab.rec.dc.cluster.FittedCluster;
import org.jlab.rec.dc.segment.Segment;
import org.jlab.rec.dc.cross.Cross;
import org.jlab.rec.dc.cross.CrossList;
import org.jlab.rec.dc.cross.CrossListFinder;
import org.jlab.rec.dc.track.Track;
import org.jlab.rec.dc.track.TrackCandListFinder;

/**
 * An engine to run the Kalman filtering process for the DC software.
 * @author benkel
 */
public class DCKFEngine extends ReconstructionEngine {

    DCGeant4Factory dcDetector;
    private int eventCounter = 0;
    private boolean debug = true;

    public DCKFEngine() {
        super("DCKF", "benkel", "0.11");
    }

    @Override
    public boolean init() {
        Constants.Load();

        // Load tables
        String[] dcTables = new String[]{
            // "/calibration/dc/signal_generation/doca_resolution",
            // "/calibration/dc/time_to_distance/t2d",
            "/calibration/dc/time_to_distance/time2dist",
            // "/calibration/dc/time_corrections/T0_correction",
            // "/calibration/dc/time_corrections/tdctimingcuts",
            // "/calibration/dc/time_jitter",
            // "/calibration/dc/tracking/wire_status",
        };

        // Load constants
        requireConstants(Arrays.asList(dcTables));

        // Load geometry
        String geomDBVar = this.getEngineConfigString("geomDBVariation");
        ConstantProvider provider =
                GeometryFactory.getConstants(DetectorType.DC,
                                             11,
                                             Optional.ofNullable(geomDBVar).orElse("default"));
        dcDetector = new DCGeant4Factory(provider, DCGeant4Factory.MINISTAGGERON);

        return true;
    }

    @Override
    public boolean processDataEvent(DataEvent event) {
        int currentEvent = eventCounter;
        eventCounter++;

        // if (currentEvent != 136) return true;

        // === INITIAL CHECKUP =========================================================
        if (!event.hasBank("RUN::config")) return true;
        DataBank headerBank = event.getBank("RUN::config");

        int newRun = headerBank.getInt("run", 0);
        if (newRun == 0) return true;

        Swim dcSwim = new Swim();
        RecoBankReader rbr = new RecoBankReader();
        RecoBankWriter rbw = new RecoBankWriter();

        // === GET OBJECTS =============================================================
        if (!event.hasBank("HitBasedTrkg::HBCrosses")) return true;

        DataBank hiBank = event.getBank("HitBasedTrkg::HBHits");
        DataBank clBank = event.getBank("HitBasedTrkg::HBClusters");
        DataBank seBank = event.getBank("HitBasedTrkg::HBSegments");
        DataBank crBank = event.getBank("HitBasedTrkg::HBCrosses");

        // Pull the hits, clusters, segments and crosses from the bank.
        if (hiBank.rows() == 0) return true;
        if (clBank.rows() == 0) return true;
        if (seBank.rows() == 0) return true;
        if (crBank.rows() == 0) return true;

        List<FittedHit>     hits     = new ArrayList();
        List<FittedCluster> clusters = new ArrayList();
        List<Segment>       segments = new ArrayList();
        List<Cross>         crosses  = new ArrayList();

        for (int h = 0; h < hiBank.rows(); ++h) hits.add    (rbr.getHit    (hiBank, h));
        for (int c = 0; c < clBank.rows(); ++c) clusters.add(rbr.getCluster(clBank, hits, c));
        for (int s = 0; s < seBank.rows(); ++s) segments.add(rbr.getSegment(seBank, clusters, s));
        for (int c = 0; c < crBank.rows(); ++c) crosses.add (rbr.getCross  (crBank, segments, c));

        // TODO: after all is up and running I should check what data from the
        //       crosses/segments/clusters/hits I can ignore so that I minimize
        //       what RecoBankReader has to read and accelerate the whole ordeal

        // === CREATE CROSSLIST FROM CROSSES ===========================================
        CrossListFinder crossLister = new CrossListFinder();

        CrossList crosslist = crossLister.candCrossLists(crosses,
                false,
                this.getConstantsManager().getConstants(newRun, Constants.TIME2DIST),
                dcDetector,
                null,
                dcSwim);

        // === INSTANCE AND RUN TrackCandListFinder ====================================
        TrackCandListFinder trkCandFinder = new TrackCandListFinder(Constants.HITBASE);
        List<Track> trkcands = trkCandFinder.getTrackCands(crosslist,
                                                           dcDetector,
                                                           Swimmer.getTorScale(),
                                                           dcSwim);

        // === WRITE TO THE BANKS ======================================================
        if (trkcands.size() > 0) {
            // TODO: Depending on how much this method takes in relation to how much writing the
            //       tracks to the bank takes, it might be more efficient to just write everything
            //       and to remove overlapping tracks in DCHB2.
            trkCandFinder.removeOverlappingTracks(trkcands);
        }

        rbw.fillAllHBBanks(event, rbw, hits, clusters, segments, crosses, trkcands);
        // rbw.fillTrajectoryBank(event, rbw, trkcands);

        return true;

    }
}
