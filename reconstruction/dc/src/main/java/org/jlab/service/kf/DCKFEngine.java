package org.jlab.service.kf;

import java.util.Arrays;
import java.util.ArrayList;
import java.util.List;
import java.util.Optional;

import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.geom.prim.Point3D;
import org.jlab.geom.base.ConstantProvider;
import org.jlab.detector.base.DetectorType;
import org.jlab.detector.base.GeometryFactory;
import org.jlab.detector.geant4.v2.DCGeant4Factory;
import org.jlab.clas.swimtools.Swim;
import org.jlab.clas.swimtools.Swimmer;
import org.jlab.clas.reco.ReconstructionEngine;
import org.jlab.rec.dc.Constants;
import org.jlab.rec.dc.banks.RecoBankWriter;
import org.jlab.rec.dc.banks.RecoBankReader;
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
    private int eventCounter;
    private int debug = 2;

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

        eventCounter = -1;
        return true;
    }

    @Override
    public boolean processDataEvent(DataEvent event) {
        eventCounter++;
        if (debug >= 1) {
            System.out.println("[DCKF] Summary for DataEvent number " +
                               eventCounter + ":");
        }

// === INITIAL CHECKUP =========================================================
        if (!event.hasBank("RUN::config")) return true;
        DataBank headerBank = event.getBank("RUN::config");

        int newRun = headerBank.getInt("run", 0);
        if (newRun == 0) return true;

        Swim dcSwim = new Swim();
        RecoBankWriter rbw = new RecoBankWriter();
        RecoBankReader rbr = new RecoBankReader();
        if (debug >= 1) {
            System.out.println("[DCKF] " + eventCounter + ".1: initialization done.");
        }
// === GET CROSSES =============================================================
        if (!event.hasBank("HitBasedTrkg::HBCrosses")  ||
            !event.hasBank("HitBasedTrkg::HBSegments") ||
            !event.hasBank("HitBasedTrkg::HBClusters") ||
            !event.hasBank("HitBasedTrkg::HBHits")) {

            System.out.println("[DCKF] ERROR: The banks were not written " +
                               "correctly in the DCHB1 engine.");
            return true;
        }
        DataBank crossesBank  = event.getBank("HitBasedTrkg::HBCrosses");
        DataBank segmentsBank = event.getBank("HitBasedTrkg::HBSegments");
        DataBank clustersBank = event.getBank("HitBasedTrkg::HBClusters");
        DataBank hitsBank     = event.getBank("HitBasedTrkg::HBHits");
        if (eventCounter == 1 && debug >= 2) {
            crossesBank.show();
            segmentsBank.show();
            clustersBank.show();
            hitsBank.show();
        }

        if (debug >= 1) {
            System.out.println("[DCKF] " + eventCounter + ".2: banks loaded.");
        }
        // Pull the crosses from the bank.
        int nCrosses = crossesBank.rows();
        if (nCrosses == 0) {
            System.out.println("[DCKF] ERROR: No crosses found in " +
                               "HitBasedTrkg::HBCrosses bank.");
            return true;
        }
        List<Cross> crosses = new ArrayList();

        // TODO: after all is up and running I should check what data from the
        //       crosses/segments/clusters/hits I can ignore so that I minimize
        //       what RecoBankReader has to read and accelerate the whole ordeal
        for (int c = 0; c < nCrosses; c++) {
            crosses.add(rbr.getCross(crossesBank, segmentsBank,
                                     clustersBank, hitsBank, c));

            if (eventCounter == 1 && debug >= 2) {
                /**************************************************************/
                System.out.println(
                    "      Cross " + crosses.get(c).get_Id() + " retrieved. Data:\n" +
                    "        " + crosses.get(c).printInfo() + "\n\n" +
                    /**********************************************************/
                    "      Segments " + crosses.get(c).get_Segment1().get_Id() +
                    " and "           + crosses.get(c).get_Segment2().get_Id() + " retrieved. Data:\n" +
                    "         " + crosses.get(c).get_Segment1().printInfo() + "\n" +
                    "         " + crosses.get(c).get_Segment2().printInfo() + "\n\n" +
                    /**********************************************************/
                    "      Clusters " + crosses.get(c).get_Segment1().get_fittedCluster().get_Id() +
                    " and "           + crosses.get(c).get_Segment2().get_fittedCluster().get_Id() +
                    " retrieved. Data:\n" +
                    "         " + crosses.get(c).get_Segment1().get_fittedCluster().printInfo() + "\n" +
                    "         " + crosses.get(c).get_Segment2().get_fittedCluster().printInfo() + "\n"
                );
                /**************************************************************/
                FittedCluster cluster1 = crosses.get(c).get_Segment1().get_fittedCluster();
                FittedCluster cluster2 = crosses.get(c).get_Segment2().get_fittedCluster();
                System.out.println("       Hits in Cluster " + cluster1.get_Id() + ":");
                for (int i = 0; i < cluster1.size(); i++) {
                    System.out.println("         " + cluster1.get(i).printInfo());
                }
                System.out.println("");
                System.out.println("       Hits in Cluster " + cluster2.get_Id() + ":");
                for (int i = 0; i < cluster2.size(); i++) {
                    System.out.println("         " + cluster2.get(i).printInfo());
                }
                System.out.println("");
                /**************************************************************/
            }
        }
        if (debug >= 1) {
            System.out.println("[DCKF] " + eventCounter + ".3: crosses loaded.");
        }
        // TODO: remove crosses, segments, clusters and hits from banks.

// === CREATE CROSSLIST FROM CROSSES ===========================================
        CrossListFinder crossLister = new CrossListFinder();

        CrossList crosslist = crossLister.candCrossLists(crosses,
                false,
                this.getConstantsManager().getConstants(newRun, Constants.TIME2DIST),
                dcDetector,
                null,
                dcSwim);
        if (debug >= 1) {
            System.out.println("[DCKF] " + eventCounter + ".4: cross lists found.");
        }

// === INSTANCE AND RUN TrackCandListFinder ====================================
        TrackCandListFinder trkCandFinder = new TrackCandListFinder(Constants.HITBASE);
        List<Track> trkcands = trkCandFinder.getTrackCands(crosslist,
                                                           dcDetector,
                                                           Swimmer.getTorScale(),
                                                           dcSwim);
        if (debug >= 1) {
            System.out.println("[DCKF] " + eventCounter + ".5: tracks found.");
        }
        if (trkcands.size() > 0) {
            trkCandFinder.removeOverlappingTracks(trkcands);
            rbw.fillHBTracksBank(event, trkcands);
            rbw.fillTrackCovMatBank(event, trkcands);
            if (debug >= 1) {
                System.out.println("[DCKF] " + eventCounter + ".6: tracks " +
                                   "written on bank.");
            }
        }
        else if (debug >= 1) {
            System.out.println("[DCKF] " + eventCounter + ".6: no tracks" +
                               " found.");
        }

        if (debug >= 1) {
            System.out.println("[DCKF].7 processed evio event " +
                               eventCounter + " successfully!\n");
        }
        return true;
    }
}
