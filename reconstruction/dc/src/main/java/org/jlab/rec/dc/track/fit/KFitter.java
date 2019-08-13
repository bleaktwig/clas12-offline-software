package org.jlab.rec.dc.track.fit;

import java.util.ArrayList;
import java.util.List;

import Jama.Matrix;
import Jama.SingularValueDecomposition;
import Jama.CholeskyDecomposition;

import org.jlab.clas.swimtools.Swim;
import org.jlab.detector.geant4.v2.DCGeant4Factory;
import org.jlab.rec.dc.track.BFieldInterpolator;
import org.jlab.rec.dc.track.Track;
import org.jlab.rec.dc.track.fit.StateVecs.CovMat;
import org.jlab.rec.dc.track.fit.StateVecs.StateVec;
import org.jlab.rec.dc.track.fit.MeasVecs.MeasVec;

/**
 * @author ziegler
 * @since 08.08.2018 modified by gurjyan
 */
public class KFitter {
    public boolean setFitFailed = false;

    double[] first_sv;
    double[] diff;

    private StateVecs sv;
    private MeasVecs mv = new MeasVecs();
    public StateVec finalStateVec;
    public CovMat finalCovMat;
    public List<org.jlab.rec.dc.trajectory.StateVec> kfStateVecsAlongTrajectory;

    public int totNumIter = 30;
    private int interpolatedNumIter; // number of iterations using the interpolated BField.
    private Swim dcSwim;
    private BFieldInterpolator bField[];
    private boolean bMode; // boolean deciding whether the interpolated or real BField should be used.

    private double newChisq = Double.POSITIVE_INFINITY;
    public double chi2    = 0;
    private double chi2kf = 0;
    public int NDF        = 0;
    public int ConvStatus = 1;

    public KFitter(Track trk, DCGeant4Factory DcDetector, boolean TimeBasedUsingHBtrack,
            Swim swimmer, BFieldInterpolator[] bField) {

        this.sv     = new StateVecs(swimmer);
        this.dcSwim = swimmer;
        this.bField = bField;
        if (TimeBasedUsingHBtrack) this.initFromHB(trk, DcDetector);
        else                       this.init(trk, DcDetector);
    }

    private void initFromHB(Track trk, DCGeant4Factory DcDetector) {
        interpolatedNumIter = -1;

        mv.setMeasVecsFromHB(trk, DcDetector);
        sv.Z = new double[mv.measurements.size()];

        for (int i = 0; i < mv.measurements.size(); i++) sv.Z[i] = mv.measurements.get(i).z;
        sv.initFromHB(trk, sv.Z[0], this);
    }

    public void init(Track trk, DCGeant4Factory DcDetector) {
        interpolatedNumIter = 20;

        mv.setMeasVecs(trk, DcDetector);
        sv.Z = new double[mv.measurements.size()];

        for (int i = 0; i < mv.measurements.size(); i++) sv.Z[i] = mv.measurements.get(i).z;
        sv.init(trk, sv.Z[0], this);
    }

    public void runFitter(int sector) {
        // System.out.printf("SECTOR: %d\n", sector);
        // long start = System.nanoTime();

        this.chi2 = 0;
        this.NDF = mv.ndf;
        int svzLength = sv.Z.length;

        // boolean earlyStop = false;

        // Print initial state vector information:
        // first_sv = new double[] {sv.trackTraj.get(0).z,  sv.trackTraj.get(0).x,  sv.trackTraj.get(0).y,
        //                          sv.trackTraj.get(0).tx, sv.trackTraj.get(0).ty, sv.trackTraj.get(0).Q,
        //                          sv.trackTraj.get(0).B,
        //                          sv.trackCov.get(0).covMat.get(0,0),
        //                          sv.trackCov.get(0).covMat.get(1,1),
        //                          sv.trackCov.get(0).covMat.get(2,2),
        //                          sv.trackCov.get(0).covMat.get(3,3),
        //                          sv.trackCov.get(0).covMat.get(4,4)
        //                         };

        // Print initial covariance matrix:
        // sv.trackCov.get(0).print();

        // Print measurements:
        // for (MeasVec m : mv.measurements) m.print();

        // System.out.printf("TRACK HERE\n");
        for (int i = 1; i <= totNumIter; i++) {
            // if (sv.trackTraj.get(0) != null) {
            //     System.out.printf("%20.12f %20.12f %20.12f %20.12f %20.12f\n",
            //             sv.trackTraj.get(0).x,  sv.trackTraj.get(0).y, sv.trackTraj.get(0).tx,
            //             sv.trackTraj.get(0).ty, sv.trackTraj.get(0).Q, newChisq);
            // }
            // System.out.printf("iter %2d/%2d\n", i, totNumIter);
            // if (i == totNumIter) {
            //     sv.trackCov.get(0).print();
            //     if (sv.trackTraj.get(0) == null) return;
            //     diff = new double[]
            //             {sv.trackTraj.get(0).z - first_sv[0], sv.trackTraj.get(0).x - first_sv[1],
            //              sv.trackTraj.get(0).y - first_sv[2], sv.trackTraj.get(0).tx - first_sv[3],
            //              sv.trackTraj.get(0).ty - first_sv[4], sv.trackTraj.get(0).Q - first_sv[5],
            //              sv.trackTraj.get(0).B - first_sv[6],
            //              sv.trackCov.get(0).covMat.get(0,0) - first_sv[7],
            //              sv.trackCov.get(0).covMat.get(1,1) - first_sv[8],
            //              sv.trackCov.get(0).covMat.get(2,2) - first_sv[9],
            //              sv.trackCov.get(0).covMat.get(3,3) - first_sv[10],
            //              sv.trackCov.get(0).covMat.get(4,4) - first_sv[11],
            //             };
            // }
            this.chi2kf = 0;
            if (i <= interpolatedNumIter) bMode = true;
            else                          bMode = false;
            if (i > 1) // Get new initial state propagating back from the last state
                if (sv.transport(bMode, sector, bField, svzLength-1, 0,
                        sv.trackTraj.get(svzLength-1), sv.trackCov.get(svzLength-1))) return;
            for (int k = 0; k < svzLength - 1; k++) {
                if (sv.transport(bMode, sector, bField, k, k + 1, sv.trackTraj.get(k), sv.trackCov.get(k))) return;
                if (this.filter(k + 1)) return;
            }

            // System.out.printf("%10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f (%10.6f)\n",
            //         sv.trackTraj.get(svzLength-1).z,  sv.trackTraj.get(svzLength-1).x,  sv.trackTraj.get(svzLength-1).y,
            //         sv.trackTraj.get(svzLength-1).tx, sv.trackTraj.get(svzLength-1).ty, sv.trackTraj.get(svzLength-1).Q,
            //         sv.trackTraj.get(svzLength-1).B,  sv.trackTraj.get(svzLength-1).deltaPath, this.chi2kf);

            // sv.trackCov.get(svzLength-1).print();
            // System.out.printf("%16.12e\n", this.chi2kf);
            if (i <= 1) continue;

            if (this.chi2kf >= newChisq) {
                // double chi2diff = -(this.chi2kf-newChisq)/newChisq;
                // if (!Double.isNaN(chi2diff)) System.out.printf("  %12.5f\n", chi2diff);

                this.ConvStatus = 1;
                continue;
            }

            if (this.finalStateVec != null) {
                if (Math.abs(sv.trackTraj.get(svzLength-1).Q -this.finalStateVec.Q) < 5.e-4
                        && Math.abs(sv.trackTraj.get(svzLength-1).x -this.finalStateVec.x)  < 1.e-4
                        && Math.abs(sv.trackTraj.get(svzLength-1).y -this.finalStateVec.y)  < 1.e-4
                        && Math.abs(sv.trackTraj.get(svzLength-1).tx-this.finalStateVec.tx) < 1.e-6
                        && Math.abs(sv.trackTraj.get(svzLength-1).ty-this.finalStateVec.ty) < 1.e-6) {
                    // System.out.printf("🍺 - Track stopped early!\n");
                    // earlyStop = true;
                    i = totNumIter;
                }
            }
            this.finalStateVec = sv.trackTraj.get(svzLength - 1);
            this.finalCovMat   = sv.trackCov.get(svzLength - 1);

            // System.out.printf("iter %2d/%2d:\n", i, totNumIter);
            // double chi2diff = -(this.chi2kf-newChisq)/newChisq;
            // if (!Double.isNaN(chi2diff)) System.out.printf("  %12.5f\n", chi2diff);

            newChisq = this.chi2kf;
        }

        if (totNumIter == 1) {
            this.finalStateVec = sv.trackTraj.get(svzLength - 1);
            this.finalCovMat   = sv.trackCov.get(svzLength - 1);
        }

        // System.out.printf("state vectors:\n");
        // for (StateVec s : sv.trackTraj.values()) s.print();
        // System.out.printf("\ncovariance matrices:\n");
        // for (CovMat   S : sv.trackCov.values())  S.print();
        // System.out.printf("\nmeasurement vectors:\n");
        // for (MeasVec  m : mv.measurements)       m.print();

        this.calcFinalChisq(sector);

        // if (diff == null) return;
        // for (double diff_i : diff) System.out.printf("%20.12f ", diff_i);
        // System.out.printf("%20.12f\n", chi2);

        // if (earlyStop) System.out.printf("  chi2: %16.12f\n", chi2);

        // long end = System.nanoTime();
        // System.out.printf("%16d\n", end - start);

        // System.out.printf("%16.12f\n", chi2);
        // System.out.printf("%5d %5d\n", sv.stepSizeDuplications, sv.totalRuns);
    }

    private boolean filter(int k) {
        if (sv.trackTraj.get(k) == null
                || sv.trackCov.get(k).covMat == null
                || k >= sv.Z.length) {

            return true;
        }

        // NOTE: ==- FIRST DETERMINANT -============================================================
        // double det1 = Math.abs(getDetDirect(sv.trackCov.get(k).covMat));
        // double det1 = Math.abs(getDetSuChang(sv.trackCov.get(k).covMat));
        // double det1 = Math.abs(getDetChol(sv.trackCov.get(k).covMat));
        // double det1 = Math.abs(getDetHDirect(sv.trackCov.get(k).covMat));
        // double det1 = Math.abs(getDetHSuChang(sv.trackCov.get(k).covMat));
        double det1 = Math.abs(getDetHChol(sv.trackCov.get(k).covMat));

        // if (det1 <= 1.e-30) {
        if (det1 <= 1.e-30 || det1 >= 1.e90 || Double.isNaN(det1)) {
            // System.out.printf("track dropped (1)\n");
            return true;
        }
        // =========================================================================================

        double[] K = new double[5];
        double V   = Math.abs(mv.measurements.get(k).unc);
        double[] H = mv.H(sv.trackTraj.get(k).y,
                     mv.measurements.get(k).tilt,
                     mv.measurements.get(k).wireMaxSag,
                     mv.measurements.get(k).wireLen);

        // ==- NOTE: ALPHA DETERMINANT =============================================================
        double detAlpha = getDetAlpha(det1, V, sv.trackCov.get(k).covMat.get(0,0),
                sv.trackCov.get(k).covMat.get(0,1), sv.trackCov.get(k).covMat.get(1,1), H[1]);
        // if (detAlpha <= 1.e-30) {
        if (detAlpha <= 1.e-30 || detAlpha >= 1.e90 || Double.isNaN(detAlpha)) {
            // System.out.printf("track dropped (2)\n");
            return true;
        }
        // =========================================================================================

        // Sherman-Morrison formula applied to the filter. This is possible because the HTGH matrix
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

        // ==- NOTE: SECOND DETERMINANT -===========================================================
        // double det2 = Math.abs(getDetDirect(sv.trackCov.get(k).covMat));
        // double det2 = Math.abs(getDetSuChang(sv.trackCov.get(k).covMat));
        // double det2 = Math.abs(getDetChol(sv.trackCov.get(k).covMat));
        // double det2 = Math.abs(getDetHDirect(sv.trackCov.get(k).covMat));
        // double det2 = Math.abs(getDetHSuChang(sv.trackCov.get(k).covMat));
        // double det2 = Math.abs(getDetHChol(sv.trackCov.get(k).covMat));

        // if (det2 <= 1.e-30) {
        // if (det2 <= 1.e-30 || det2 >= 1.e60 || Double.isNaN(det2)) {
        //     System.out.printf("track dropped (2)\n");
        //     return true;
        // }
        // =========================================================================================

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

        return false;
    }

    private void calcFinalChisq(int sector) {
        int k       = sv.Z.length - 1;
        this.chi2   = 0;
        double path = 0;
        kfStateVecsAlongTrajectory = new ArrayList<>();

        if (sv.trackTraj.get(k) == null || sv.trackCov.get(k).covMat == null)
            return;

        sv.transport(false, sector, null, sv.Z.length - 1, 0,
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
            sv.transport(false, sector, null, k1 - 1, k1, sv.trackTraj.get(k1 - 1), sv.trackCov.get(k1 - 1));

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

    private double getDetDirect(Matrix mat) {
        return mat.det();
    }

    private double getDetSuChang(Matrix A) {
        int rd = A.getRowDimension();
        int cd = A.getColumnDimension();
        if (rd == 1) return A.get(0,0);

        Matrix W = A.getMatrix(0, rd-2, 0, cd-2);
        Matrix v = A.getMatrix(0, rd-2, cd-1, cd-1);
        Matrix u = A.getMatrix(rd-1, rd-1, 0, cd-2);
        double r = A.get(rd-1, cd-1);

        Matrix vu_r = v.times(u);
        for (int i = 0; i < vu_r.getRowDimension(); ++i) {
            for (int j = 0; j < vu_r.getColumnDimension(); ++j) {
                vu_r.set(i,j, vu_r.get(i,j)/r);
            }
        }

        return r * getDetSuChang(W.minus(vu_r));
    }

    private double getDetChol(Matrix mat) {
        CholeskyDecomposition c = new CholeskyDecomposition(mat);
        double diag = c.getL().get(0,0) * c.getL().get(1,1) * c.getL().get(2,2) * c.getL().get(3,3) * c.getL().get(4,4);
        return diag*diag;
    }

    private double getDetHDirect(Matrix mat) {
        double det;
        {
            double p0 = mat.get(0,1)*mat.get(2,3);
            double p1 = mat.get(0,1)*mat.get(2,4);
            double p2 = mat.get(1,3)*mat.get(2,4);
            double p3 = mat.get(1,4)*mat.get(2,3);
            double p4 = mat.get(1,2)*mat.get(3,4);
            double p5 = mat.get(0,0)*mat.get(1,1);
            double p6 = mat.get(0,2)*mat.get(3,3);
            double p7 = mat.get(0,3)*mat.get(0,3);
            double p8 = mat.get(0,4)*mat.get(0,4);
            double p9 = mat.get(1,2)*mat.get(4,4);

            double r0 = mat.get(1,3)*mat.get(2,2);
            double r1 = mat.get(2,2)*mat.get(3,3);
            double r2 = mat.get(0,2)*mat.get(1,3);
            double r3 = mat.get(0,0)*mat.get(1,4);
            double r4 = mat.get(0,3)*mat.get(0,4);
            double r5 = mat.get(1,1)*mat.get(2,3);
            double r6 = mat.get(1,2)*mat.get(3,3);
            double r7 = mat.get(0,2)*mat.get(1,4);

            double q00 = mat.get(0,0)*mat.get(1,2)*mat.get(3,4);
            double q01 = mat.get(0,1)*mat.get(0,2)*mat.get(3,4);
            double q02 = mat.get(0,2)*mat.get(0,3)*mat.get(1,4);
            double q03 = mat.get(0,3)*mat.get(0,4)*mat.get(1,2);
            double q04 = p0*mat.get(0,4);
            double q05 = p0*mat.get(4,4);
            double q06 = p1*mat.get(0,3);
            double q07 = p1*mat.get(3,3);
            double q08 = mat.get(0,1)*mat.get(2,2)*mat.get(3,4);
            double q09 = mat.get(0,2)*mat.get(0,4)*mat.get(1,3);
            double q10 = mat.get(0,2)*mat.get(1,1)*mat.get(3,4);

            det = p5*r1*mat.get(4,4)  -   p5*mat.get(2,2)*mat.get(3,4)*mat.get(3,4)
                    -   p5*mat.get(2,3)*mat.get(2,3)*mat.get(4,4) + 2*p5*mat.get(2,3)*mat.get(2,4)*mat.get(3,4)
                    -   p5*mat.get(2,4)*mat.get(2,4)*mat.get(3,3) -   p9*r6*mat.get(0,0)
                    +   q00*p4   + 2*p9*mat.get(0,0)*mat.get(1,3)*mat.get(2,3)
                    - 2*q00*p2   - 2*q00*p3
                    + 2*r3*r6*mat.get(2,4)  -   r0*mat.get(0,0)*mat.get(1,3)*mat.get(4,4)
                    +   mat.get(0,0)*p2*p2  + 2*r0*r3*mat.get(3,4)
                    - 2*mat.get(0,0)*p3*p2  -   r1*r3*mat.get(1,4)
                    +   mat.get(0,0)*p3*p3  -   r1*mat.get(0,1)*mat.get(0,1)*mat.get(4,4)
                    +   q08*mat.get(0,1)*mat.get(3,4)  +   p0*q05
                    - 2*p0*p1*mat.get(3,4)  +   p1*q07
                    + 2*p6*p9*mat.get(0,1)  - 2*q01*p4
                    - 2*q05*r2   + 2*q01*p2
                    + 2*q01*p3   - 2*q07*r7
                    - 2*mat.get(0,3)*mat.get(1,2)*q05  + 2*q06*p4
                    + 2*r0*mat.get(0,1)*mat.get(0,3)*mat.get(4,4) - 2*q06*p2
                    - 2*q08*mat.get(0,3)*mat.get(1,4)  + 2*q06*p3
                    + 2*p4*q04   - 2*q07*mat.get(0,4)*mat.get(1,2)
                    - 2*q08*mat.get(0,4)*mat.get(1,3)  + 2*q04*p2
                    + 2*r1*mat.get(0,1)*mat.get(0,4)*mat.get(1,4) - 2*q04*p3
                    -   p6*mat.get(0,2)*mat.get(1,1)*mat.get(4,4) +   q10*mat.get(0,2)*mat.get(3,4)
                    +   r2*r2*mat.get(4,4)  - 2*r2*r7*mat.get(3,4)
                    +   p6*r7*mat.get(1,4)  + 2*r5*mat.get(0,2)*mat.get(0,3)*mat.get(4,4)
                    - 2*q10*mat.get(0,3)*mat.get(2,4)  - 2*r2*mat.get(0,3)*p9
                    + 2*q02*p4   + 2*q02*p2
                    - 2*q02*p3   - 2*q10*mat.get(0,4)*mat.get(2,3)
                    + 2*p6*mat.get(0,4)*mat.get(1,1)*mat.get(2,4) + 2*q09*p4
                    - 2*p6*mat.get(0,4)*mat.get(1,2)*mat.get(1,4) - 2*q09*p2
                    + 2*q09*p3   -   p7*mat.get(1,1)*mat.get(2,2)*mat.get(4,4)
                    +   p7*mat.get(1,1)*mat.get(2,4)*mat.get(2,4) +   p7*p9*mat.get(1,2)
                    - 2*p7*mat.get(1,2)*mat.get(1,4)*mat.get(2,4) +   p7*mat.get(1,4)*mat.get(1,4)*mat.get(2,2)
                    + 2*r4*mat.get(1,1)*mat.get(2,2)*mat.get(3,4) - 2*r4*r5*mat.get(2,4)
                    - 2*q03*p4   + 2*q03*p2
                    + 2*q03*p3   - 2*r0*r4*mat.get(1,4)
                    -   p8*mat.get(1,1)*r1  +   p8*r5*mat.get(2,3)
                    +   p8*r6*mat.get(1,2)  - 2*p8*mat.get(1,2)*mat.get(1,3)*mat.get(2,3)
                    +   p8*r0*mat.get(1,3);
        }
        return det;
    }

    private double getDetHSuChang(Matrix mat) {
        double det;
        {
            double p0 = mat.get(0,0)*mat.get(4,4) - mat.get(0,4)*mat.get(0,4);
            double p1 = mat.get(0,1)*mat.get(4,4) - mat.get(0,4)*mat.get(1,4);
            double p2 = mat.get(0,2)*mat.get(4,4) - mat.get(0,4)*mat.get(2,4);
            double p3 = mat.get(1,1)*mat.get(4,4) - mat.get(1,4)*mat.get(1,4);
            double p4 = mat.get(1,2)*mat.get(4,4) - mat.get(1,4)*mat.get(2,4);
            double p5 = mat.get(2,2)*mat.get(4,4) - mat.get(2,4)*mat.get(2,4);

            double q0 = mat.get(0,3)*mat.get(4,4) - mat.get(0,4)*mat.get(3,4);
            double q1 = mat.get(1,3)*mat.get(4,4) - mat.get(1,4)*mat.get(3,4);
            double q2 = mat.get(2,3)*mat.get(4,4) - mat.get(2,4)*mat.get(3,4);
            double q3 = mat.get(3,3)*mat.get(4,4) - mat.get(3,4)*mat.get(3,4);

            double r0 = p5*q3 - q2*q2;
            double r1 = p5*q3 - q0*q2;
            double r2 = p4*q3 - q1*q2;

            det = 1/(mat.get(4,4)*mat.get(4,4)*mat.get(4,4)*q3*q3*r0)
                    * (((p0*q3 - q0*q0)*r0 - r1*r1) * ((p3*q3 - q1*q1)*r0 - r2) - ((p1*q3 - q0*q1)*r0 - r1*r2));
        }
        return det;
    }

    private double getDetHChol(Matrix mat) {
        double q00 = Math.sqrt(mat.get(0,0));
        double q10 = mat.get(0,1)/q00;
        double q20 = mat.get(0,2)/q00;
        double q30 = mat.get(0,3)/q00;
        double q40 = mat.get(0,4)/q00;

        double q11 = Math.sqrt(mat.get(1,1) - q10*q10);
        double q21 = (mat.get(1,2) - q10*q20) / q11;
        double q31 = (mat.get(1,3) - q10*q30) / q11;
        double q41 = (mat.get(1,4) - q10*q40) / q11;

        double q22 = Math.sqrt(mat.get(2,2) - q21*q21 - q20*q20);
        double q32 = (mat.get(2,3) - q21*q31 - q20*q30) / q22;
        double q42 = (mat.get(2,4) - q21*q41 - q20*q40) / q22;

        double q33 = Math.sqrt(mat.get(3,3) - q32*q32 - q31*q31 - q30*q30);
        double q43 = (mat.get(3,4) - q32*q42 - q31*q41 - q30*q40) / q33;

        double q44 = Math.sqrt(mat.get(4,4) - q43*q43 - q42*q42 - q41*q41 - q40*q40);

        double sqrtDet = q00 * q11 * q22 * q33 * q44;
        return sqrtDet*sqrtDet;
    }

    private double getDetAlpha(double detCk_1, double v, double a, double b, double f, double p) {
        return (detCk_1*v)/(f*p*p + 2*b*p + a + v);
    }
}
