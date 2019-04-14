package org.jlab.rec.dc.track.fit;

import java.util.ArrayList;
import java.util.List;

import Jama.Matrix;

import org.jlab.clas.swimtools.Swim;
import org.jlab.detector.geant4.v2.DCGeant4Factory;
import org.jlab.rec.dc.track.Track;
import org.jlab.rec.dc.track.fit.StateVecs.CovMat;
import org.jlab.rec.dc.track.fit.StateVecs.StateVec;


/**
 * @author ziegler
 * @since 08.08.2018 modified by gurjyan
 */
public class KFitter {
    public boolean setFitFailed = false;

    private StateVecs sv;
    private MeasVecs mv = new MeasVecs();

    public StateVec finalStateVec;
    public CovMat finalCovMat;
    public List<org.jlab.rec.dc.trajectory.StateVec> kfStateVecsAlongTrajectory;
    public int totNumIter   = 30;
    private double newChisq = Double.POSITIVE_INFINITY;

    public double chi2    = 0;
    private double chi2kf = 0;
    public int NDF        = 0;
    public int ConvStatus = 1;

    public KFitter(Track trk, DCGeant4Factory DcDetector, boolean TimeBasedUsingHBtrack,
            Swim swimmer) {

        sv = new StateVecs(swimmer);
        if (TimeBasedUsingHBtrack) this.initFromHB(trk, DcDetector);
        else                       this.init(trk, DcDetector);
    }

    private void initFromHB(Track trk, DCGeant4Factory DcDetector) {
        mv.setMeasVecsFromHB(trk, DcDetector);
        sv.Z = new double[mv.measurements.size()];

        for (int i = 0; i < mv.measurements.size(); i++) sv.Z[i] = mv.measurements.get(i).z;
        sv.initFromHB(trk, sv.Z[0], this);
    }

    public void init(Track trk, DCGeant4Factory DcDetector) {
        mv.setMeasVecs(trk, DcDetector);
        sv.Z = new double[mv.measurements.size()];

        for (int i = 0; i < mv.measurements.size(); i++) sv.Z[i] = mv.measurements.get(i).z;
        sv.init(trk, sv.Z[0], this);
    }

    public void runFitter(int sector) {
        this.chi2 = 0;
        this.NDF = mv.ndf;
        int svzLength = sv.Z.length;

        for (int i = 1; i <= totNumIter; i++) {
            this.chi2kf = 0;
            if (i > 1) {
                // Get new state vec at 1st measurement after propagating back from the last
                // filtered state
                sv.transport(sector, svzLength-1, 0, sv.trackTraj.get(svzLength-1),
                        sv.trackCov.get(svzLength-1));
            }
            for (int k = 0; k < svzLength - 1; k++) {
                sv.transport(sector, k, k + 1, sv.trackTraj.get(k), sv.trackCov.get(k));
                this.filter(k + 1);
            }
            if (i <= 1) continue;

            if (this.chi2kf >= newChisq) {
                this.ConvStatus = 1;
                continue;
            }

            if (this.finalStateVec != null) {
                if (Math.abs(sv.trackTraj.get(svzLength-1).Q -this.finalStateVec.Q) < 5.e-4
                        && Math.abs(sv.trackTraj.get(svzLength-1).x -this.finalStateVec.x)  < 1.e-4
                        && Math.abs(sv.trackTraj.get(svzLength-1).y -this.finalStateVec.y)  < 1.e-4
                        && Math.abs(sv.trackTraj.get(svzLength-1).tx-this.finalStateVec.tx) < 1.e-6
                        && Math.abs(sv.trackTraj.get(svzLength-1).ty-this.finalStateVec.ty) < 1.e-6) {
                    i = totNumIter;
                }
            }
            this.finalStateVec = sv.trackTraj.get(svzLength - 1);
            this.finalCovMat   = sv.trackCov.get(svzLength - 1);

            newChisq = this.chi2kf;
        }

        if (totNumIter == 1) {
            this.finalStateVec = sv.trackTraj.get(svzLength - 1);
            this.finalCovMat   = sv.trackCov.get(svzLength - 1);
        }

        this.calcFinalChisq(sector);
    }

    private void filter(int k) {
        if (sv.trackTraj.get(k)       == null ||
            sv.trackCov.get(k).covMat == null ||
            k                         >= sv.Z.length) {

            return;
        }

        double[] K = new double[5];
        double V   = Math.abs(mv.measurements.get(k).unc);
        double[] H = mv.H(sv.trackTraj.get(k).y,
                     mv.measurements.get(k).tilt,
                     mv.measurements.get(k).wireMaxSag,
                     mv.measurements.get(k).wireLen);

        // Sherman-Morrisey formula applied to the filter. This is possible because the HTGH matrix
        // used can be expressed as a column vector multiplied by a row vector (in this case, itself
        // transposed).
        Matrix Hvec  = new Matrix(new double[][] {{H[0]}, {H[1]}, {0}, {0}, {0}});
        Matrix HvecT = Hvec.transpose();
        Matrix C     = sv.trackCov.get(k).covMat;

        double div = (new Matrix(new double[][] {{V}}).plus((HvecT.times(C)).times(Hvec))).get(0, 0);
        Matrix result = ((C.times(Hvec)).times(HvecT)).times(C);
        for (int i = 0; i < 5; ++i) {
            for (int j = 0; j < 5; ++j) {
                result.set(i, j, result.get(i, j)/div);
            }
        }
        sv.trackCov.get(k).covMat = C.minus(result);

        for (int j = 0; j < 5; j++) {
            // the gain matrix
            K[j] = (H[0] * sv.trackCov.get(k).covMat.get(j, 0) +
                    H[1] * sv.trackCov.get(k).covMat.get(j, 1)) / V;
        }

        double h = mv.h(new double[]{sv.trackTraj.get(k).x, sv.trackTraj.get(k).y},
                        mv.measurements.get(k).tilt,
                        mv.measurements.get(k).wireMaxSag,
                        mv.measurements.get(k).wireLen);

        chi2kf += ((mv.measurements.get(k).x - h) * (mv.measurements.get(k).x - h) / V);
        sv.trackTraj.get(k).x  += K[0] * (mv.measurements.get(k).x - h);
        sv.trackTraj.get(k).y  += K[1] * (mv.measurements.get(k).x - h);
        sv.trackTraj.get(k).tx += K[2] * (mv.measurements.get(k).x - h);
        sv.trackTraj.get(k).ty += K[3] * (mv.measurements.get(k).x - h);
        sv.trackTraj.get(k).Q  += K[4] * (mv.measurements.get(k).x - h);
    }

    @SuppressWarnings("unused")
    private void smooth(int sector, int k) {
        this.chi2 = 0;
        if (sv.trackTraj.get(k) == null || sv.trackCov.get(k).covMat == null)
            return;

        sv.transport(sector, k, 0, sv.trackTraj.get(k), sv.trackCov.get(k));
        for (int k1 = 0; k1 < k; k1++) {
            sv.transport(sector, k1, k1 + 1, sv.trackTraj.get(k1), sv.trackCov.get(k1));
            this.filter(k1 + 1);
        }
    }

    private void calcFinalChisq(int sector) {
        int k       = sv.Z.length - 1;
        this.chi2   = 0;
        double path = 0;
        kfStateVecsAlongTrajectory = new ArrayList<>();

        if (sv.trackTraj.get(k) == null || sv.trackCov.get(k).covMat == null)
            return;

        sv.transport(sector, sv.Z.length - 1, 0,
                     sv.trackTraj.get(sv.Z.length - 1),
                     sv.trackCov.get(sv.Z.length - 1));

        org.jlab.rec.dc.trajectory.StateVec svc =
            new org.jlab.rec.dc.trajectory.StateVec(sv.trackTraj.get(0).x,
                                                    sv.trackTraj.get(0).y,
                                                    sv.trackTraj.get(0).tx,
                                                    sv.trackTraj.get(0).ty);

        svc.setZ(sv.trackTraj.get(0).z);
        svc.setB(sv.trackTraj.get(0).B);
        path += sv.trackTraj.get(0).deltaPath;
        svc.setPathLength(path);
        double h0 = mv.h(new double[]{sv.trackTraj.get(0).x, sv.trackTraj.get(0).y},
                         mv.measurements.get(0).tilt,
                         mv.measurements.get(0).wireMaxSag,
                         mv.measurements.get(0).wireLen);

        svc.setProjector(h0);
        kfStateVecsAlongTrajectory.add(svc);

        chi2 += (mv.measurements.get(0).x - h0) * (mv.measurements.get(0).x - h0) /
                mv.measurements.get(0).error;

        for (int k1 = 1; k1 <= k; k1++) {
            sv.transport(sector, k1 - 1, k1, sv.trackTraj.get(k1 - 1), sv.trackCov.get(k1 - 1));

            double V = mv.measurements.get(k1).error;
            double h = mv.h(new double[]{sv.trackTraj.get(k1).x, sv.trackTraj.get(k1).y},
                            mv.measurements.get(k1).tilt,
                            mv.measurements.get(k1).wireMaxSag,
                            mv.measurements.get(k1).wireLen);

            svc = new org.jlab.rec.dc.trajectory.StateVec(sv.trackTraj.get(k1).x,
                                                          sv.trackTraj.get(k1).y,
                                                          sv.trackTraj.get(k1).tx,
                                                          sv.trackTraj.get(k1).ty);

            svc.setZ(sv.trackTraj.get(k1).z);
            svc.setB(sv.trackTraj.get(k1).B);
            path += sv.trackTraj.get(k1).deltaPath;
            svc.setPathLength(path);
            svc.setProjector(h);
            kfStateVecsAlongTrajectory.add(svc);
            chi2 += (mv.measurements.get(k1).x - h) * (mv.measurements.get(k1).x - h) / V;
        }
    }

    private boolean isNonsingular(Matrix mat) {
        return Math.abs(mat.det()) >= 1.e-30;
    }
}
