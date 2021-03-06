package org.jlab.service.dc;

import java.util.ArrayList;
import java.util.List;

import Jama.Matrix;

import org.jlab.clas.swimtools.Swim;
import org.jlab.geom.prim.Point3D;
import org.jlab.geom.prim.Vector3D;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;

import org.jlab.rec.dc.Constants;
import org.jlab.rec.dc.banks.HitReader;
import org.jlab.rec.dc.banks.RecoBankWriter;
import org.jlab.rec.dc.cluster.ClusterCleanerUtilities;
import org.jlab.rec.dc.cluster.ClusterFinder;
import org.jlab.rec.dc.cluster.ClusterFitter;
import org.jlab.rec.dc.cluster.FittedCluster;
import org.jlab.rec.dc.cross.Cross;
import org.jlab.rec.dc.cross.CrossMaker;
import org.jlab.rec.dc.hit.FittedHit;
import org.jlab.rec.dc.segment.Segment;
import org.jlab.rec.dc.segment.SegmentFinder;
import org.jlab.rec.dc.timetodistance.TimeToDistanceEstimator;
import org.jlab.rec.dc.track.BFieldInterpolator;
import org.jlab.rec.dc.track.Track;
import org.jlab.rec.dc.track.TrackCandListFinder;
import org.jlab.rec.dc.track.fit.KFitter;
import org.jlab.rec.dc.trajectory.StateVec;
import org.jlab.rec.dc.trajectory.Trajectory;
import org.jlab.rec.dc.trajectory.TrajectoryFinder;

/**
 * The DC Time-based engine.
 * @author ziegler
 */
public class DCTBEngine extends DCEngine {

    private TimeToDistanceEstimator tde;
    private double tarCent = -1.942;

    private BFieldInterpolator[] bField;

    public DCTBEngine() {
        super("DCTB");
        tde = new TimeToDistanceEstimator();
    }

    @Override
    public boolean init() {
        super.LoadTables();
        return true;
    }

    @Override
    public boolean processDataEvent(DataEvent event) {
        if (event.hasBank("RUN::config") == false) {
            System.err.println("RUN CONDITIONS NOT READ AT TIMEBASED LEVEL!");
            return true;
        }
        if (event.hasBank("MC::Event") == true) tarCent = 0;
        // if (event.getBank("RECHB::Event").getFloat("STTime", 0) < 0)
        //     return true; // Require the start time to reconstruct the tracks in the event

        DataBank bank = event.getBank("RUN::config");

        // Load the constants
        int newRun = bank.getInt("run", 0);
        if (newRun == 0) return true;

        double T_Start = 0;
        if (Constants.isUSETSTART()) {
            if (event.hasBank("RECHB::Event")) {
                T_Start = event.getBank("RECHB::Event").getFloat("STTime", 0);
                // Quit if start time not found in data
                if (T_Start < 0) return true;
            }
            else return true; // No REC HB bank
        }

        // Get field
        Swim dcSwim                = new Swim();
        ClusterFitter cf           = new ClusterFitter();
        ClusterCleanerUtilities ct = new ClusterCleanerUtilities();

        List<FittedHit> fhits        = new ArrayList<FittedHit>();
        List<FittedCluster> clusters = new ArrayList<FittedCluster>();
        List<Segment> segments       = new ArrayList<Segment>();
        List<Cross> crosses          = new ArrayList<Cross>();
        List<Track> trkcands         = new ArrayList<Track>();

        // Instantiate bank writer
        RecoBankWriter rbc = new RecoBankWriter();

        HitReader hitRead  = new HitReader();
        hitRead.readHBHits(
                event,
                super.getConstantsManager().getConstants(newRun,
                                                "/calibration/dc/signal_generation/doca_resolution"),
                super.getConstantsManager().getConstants(newRun,
                                                "/calibration/dc/time_to_distance/time2dist"),
                Constants.getT0(), Constants.getT0Err(), dcDetector, tde);

        hitRead.readTBHits(
                event,
                super.getConstantsManager().getConstants(newRun,
                                                "/calibration/dc/signal_generation/doca_resolution"),
                super.getConstantsManager().getConstants(newRun,
                                                "/calibration/dc/time_to_distance/time2dist"),
                tde, Constants.getT0(), Constants.getT0Err()
        );

        List<FittedHit> hits = new ArrayList<FittedHit>();
        // I) Get the hits =============================================================
        hits = hitRead.get_HBHits();

        // II) Process the hits ========================================================
        // 1)  Exit if hit list is empty -----------------------------------------------
        if(hits.isEmpty()) return true;

        // 2)  Find the clusters from these hits ---------------------------------------
        ClusterFinder clusFinder = new ClusterFinder();

        clusters = clusFinder.FindTimeBasedClusters(hits, cf, ct,
            super.getConstantsManager()
                 .getConstants(newRun, "/calibration/dc/time_to_distance/time2dist"), dcDetector, tde);

        if (clusters.isEmpty()) {
            rbc.fillAllTBBanks(event, rbc, hits, null, null, null, null);
            return true;
        }

        // 3)  Find the segments from the fitted clusters ------------------------------
        SegmentFinder segFinder = new SegmentFinder();

        List<FittedCluster> pclusters = segFinder.selectTimeBasedSegments(clusters);

        segments = segFinder.get_Segments(pclusters, event, dcDetector, false);

        if (segments.isEmpty()) { // 6 segments are needed to make a trajectory
            for (FittedCluster c : clusters) {
                for (FittedHit hit : c) {
                    hit.set_AssociatedClusterID(c.get_Id());
                    hit.set_AssociatedHBTrackID(c.get(0).get_AssociatedHBTrackID());
                    fhits.add(hit);
                }
            }
            rbc.fillAllTBBanks(event, rbc, fhits, clusters, null, null, null);
            return true;
        }

        for (Segment seg : segments) {
            for (FittedHit hit : seg.get_fittedCluster()) {
                fhits.add(hit);
            }
        }

        CrossMaker crossMake = new CrossMaker();

        // The track bank is also needed
        if (!event.hasBank("HitBasedTrkg::HBTracks")) return true;

        DataBank trkbank    = event.getBank("HitBasedTrkg::HBTracks");
        DataBank trkcovbank = event.getBank("TimeBasedTrkg::TBCovMat");
        int trkrows = trkbank.rows();

        if (trkbank.rows() != trkcovbank.rows()) return true; // HB tracks not saved correctly

        Track[] TrackArray = new Track[trkrows];
        for (int i = 0; i < trkrows; i++) {
            Track HBtrk = new Track();
            HBtrk.set_Id(trkbank.getShort("id", i));
            HBtrk.set_Sector(trkbank.getByte("sector", i));
            HBtrk.set_Q(trkbank.getByte("q", i));
            HBtrk.set_pAtOrig(new Vector3D(trkbank.getFloat("p0_x", i),
                                           trkbank.getFloat("p0_y", i),
                                           trkbank.getFloat("p0_z", i)));

            HBtrk.set_Vtx0(new Point3D(trkbank.getFloat("Vtx0_x", i),
                                       trkbank.getFloat("Vtx0_y", i),
                                       trkbank.getFloat("Vtx0_z", i)));

            HBtrk.set_FitChi2(trkbank.getFloat("chi2", i));
            Matrix initCMatrix = new Matrix(new double[][]{
                {trkcovbank.getFloat("C11", i), trkcovbank.getFloat("C12", i),
                 trkcovbank.getFloat("C13", i), trkcovbank.getFloat("C14", i),
                 trkcovbank.getFloat("C15", i)},
                {trkcovbank.getFloat("C21", i), trkcovbank.getFloat("C22", i),
                 trkcovbank.getFloat("C23", i), trkcovbank.getFloat("C24", i),
                 trkcovbank.getFloat("C25", i)},
                {trkcovbank.getFloat("C31", i), trkcovbank.getFloat("C32", i),
                 trkcovbank.getFloat("C33", i), trkcovbank.getFloat("C34", i),
                 trkcovbank.getFloat("C35", i)},
                {trkcovbank.getFloat("C41", i), trkcovbank.getFloat("C42", i),
                 trkcovbank.getFloat("C43", i), trkcovbank.getFloat("C44", i),
                 trkcovbank.getFloat("C45", i)},
                {trkcovbank.getFloat("C51", i), trkcovbank.getFloat("C52", i),
                 trkcovbank.getFloat("C53", i), trkcovbank.getFloat("C54", i),
                 trkcovbank.getFloat("C55", i)}
            });

            HBtrk.set_CovMat(initCMatrix);
            TrackArray[HBtrk.get_Id()-1] = HBtrk;
            TrackArray[HBtrk.get_Id()-1].set_Status(0);
        }
        if (TrackArray == null) return true; // HB tracks not saved correctly
        for (Segment seg : segments) {
            TrackArray[seg.get(0).get_AssociatedHBTrackID()-1].get_ListOfHBSegments().add(seg);
            if (seg.get_Status() == 1)
                TrackArray[seg.get(0).get_AssociatedHBTrackID()-1].set_Status(1);
        }

        // 4)  Find the list of track candidates ---------------------------------------
        TrackCandListFinder trkcandFinder = new TrackCandListFinder("TimeBased");
        TrajectoryFinder trjFind = new TrajectoryFinder();
        for (int i = 0; i < TrackArray.length; i++) {
            if (TrackArray[i] == null
                    || TrackArray[i].get_ListOfHBSegments() == null
                    || TrackArray[i].get_ListOfHBSegments().size() < 4) {
                continue;
            }
            TrackArray[i].set_MissingSuperlayer(get_Status(TrackArray[i]));
            TrackArray[i].addAll(crossMake.find_Crosses(TrackArray[i].get_ListOfHBSegments(),
                                                        dcDetector));
            if (TrackArray[i].size() < 1) continue;
            crosses.addAll(TrackArray[i]);

            KFitter kFit = new KFitter(TrackArray[i], dcDetector, true, dcSwim, null);

            StateVec fn = new StateVec();
            kFit.runFitter(TrackArray[i].get(0).get_Sector());

            if (!kFit.setFitFailed && kFit.finalStateVec != null) {
                // Set the state vector at the last measurement site
                fn.set(kFit.finalStateVec.x,  kFit.finalStateVec.y,
                       kFit.finalStateVec.tx, kFit.finalStateVec.ty);

                // Set the track parameters if the filter does not fail
                TrackArray[i].set_P(1./Math.abs(kFit.finalStateVec.Q));
                TrackArray[i].set_Q((int)Math.signum(kFit.finalStateVec.Q));
                trkcandFinder.setTrackPars(TrackArray[i], new Trajectory(), trjFind, fn,
                        kFit.finalStateVec.z, dcDetector, dcSwim);

                // Candidate parameters are set from the state vector
                TrackArray[i].set_FitChi2(kFit.chi2);
                TrackArray[i].set_FitNDF(kFit.NDF);
                TrackArray[i].set_Trajectory(kFit.kfStateVecsAlongTrajectory);
                TrackArray[i].set_FitConvergenceStatus(kFit.ConvStatus);
                TrackArray[i].set_Id(TrackArray[i].size()+1);
                TrackArray[i].set_CovMat(kFit.finalCovMat.covMat);

                if (TrackArray[i].get_Vtx0().toVector3D().mag() > 500) continue;
                trkcands.add(TrackArray[i]);
            }
        }

        for (int i = 0; i < crosses.size(); i++) {
            crosses.get(i).set_Id(i + 1);
        }

        // Track found
        int trkId = 1;

        if (trkcands.size() > 0) {
            for (Track trk : trkcands) {
                // Reset the id
                trk.set_Id(trkId);
                trkcandFinder.matchHits(trk.get_Trajectory(),
                                        trk, dcDetector, dcSwim);
                trk.calcTrajectory(trkId, dcSwim,
                                   trk.get_Vtx0().x(),
                                   trk.get_Vtx0().y(),
                                   trk.get_Vtx0().z(),
                                   trk.get_pAtOrig().x(),
                                   trk.get_pAtOrig().y(),
                                   trk.get_pAtOrig().z(),
                                   trk.get_Q(),
                                   ftofDetector, tSurf, tarCent);

                for(Cross c : trk) {
                    c.get_Segment1().isOnTrack=true;
                    c.get_Segment2().isOnTrack=true;

                    for(FittedHit h1 : c.get_Segment1()) {
                        h1.set_AssociatedTBTrackID(trk.get_Id());
                    }
                    for (FittedHit h1 : c.get_Segment1()) {
                        h1.set_AssociatedTBTrackID(trk.get_Id());
                    }
                    for (FittedHit h2 : c.get_Segment2()) {
                        h2.set_AssociatedTBTrackID(trk.get_Id());
                    }
                }
                trkId++;
            }
        }

        if (trkcands.isEmpty()) {
            // No candidates found, stop here and save the hits, the clusters, the segments and the
            //     crosses.
            rbc.fillAllTBBanks(event, rbc, fhits, clusters, segments, crosses, null);
        }
        else {
            rbc.fillAllTBBanks(event, rbc, fhits, clusters, segments, crosses, trkcands);
        }
        return true;
    }

    private int get_Status(Track track) {
        int miss = 0;

        int L[] = new int[6];
        for(int l = 0; l < track.get_ListOfHBSegments().size(); l++) {
            L[track.get_ListOfHBSegments().get(l).get_Superlayer() - 1]++;
        }
        for(int l = 0; l < 6; l++) {
            if(L[l] == 0) miss = l + 1;
        }
        return miss;
    }
}
