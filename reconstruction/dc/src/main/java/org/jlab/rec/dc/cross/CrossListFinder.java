package org.jlab.rec.dc.cross;

import java.util.ArrayList;
import java.util.List;
import Jama.Matrix;

import org.jlab.geom.prim.Point3D;
import org.jlab.geom.prim.Vector3D;
import org.jlab.clas.swimtools.Swim;
import trackfitter.fitter.LineFitter;
import org.jlab.utils.groups.IndexedTable;
import org.jlab.detector.geant4.v2.DCGeant4Factory;

import org.jlab.rec.dc.Constants;
import org.jlab.rec.dc.hit.FittedHit;
import org.jlab.rec.dc.cluster.ClusterFitter;
import org.jlab.rec.dc.segment.Segment;
import org.jlab.rec.dc.segment.SegmentFinder;
import org.jlab.rec.dc.timetodistance.TimeToDistanceEstimator;

/**
 * A class with methods used to find lists of crosses. This is the Pattern Recognition step used in
 * track seeding, to find the points that are consistent with belonging to the same track. This step
 * precedes the initial estimates of the track parameters.
 *
 * @author ziegler
 */
public class CrossListFinder  {

    private final List<BaseCand> trkCnds = new ArrayList<BaseCand>();
    ClusterFitter cf = new ClusterFitter();
    SegmentFinder segFinder = new SegmentFinder();

    /**
     * Determines which crosses from a list are consistent with belonging to a track in the DC.
     * @param dccrosslist the list of crosses in the event
     * @param TimeBased   flag telling if the detection is hit-based or time-based
     * @param tab
     * @param DcDetector  DC detector utility
     * @param tde
     * @param swimmer
     * @return the list of crosses determined to be consistent with belonging to a track in the DC
     */
    public CrossList candCrossLists(List<Cross> dccrosslist, boolean TimeBased, IndexedTable tab,
            DCGeant4Factory DcDetector, TimeToDistanceEstimator tde, Swim swimmer) {
        trkCnds.clear();

        if (dccrosslist.size() <= 0) return new CrossList();

        List<Cross> dccrosslistRg1 = new ArrayList<Cross>();
        List<Cross> dccrosslistRg2 = new ArrayList<Cross>();
        List<Cross> dccrosslistRg3 = new ArrayList<Cross>();
        for (Cross dc : dccrosslist) {
            if (dc.get_Region() == 1) dccrosslistRg1.add(dc);
            if (dc.get_Region() == 2) dccrosslistRg2.add(dc);
            if (dc.get_Region() == 3) dccrosslistRg3.add(dc);
        }

        // Need 3 crosses
        if (!dccrosslistRg1.isEmpty()
                && !dccrosslistRg2.isEmpty()
                && !dccrosslistRg3.isEmpty()) {
            for (Cross c1 : dccrosslistRg1) {
                for (Cross c2 : dccrosslistRg2) {
                    for (Cross c3 : dccrosslistRg3) {
                        if (c1.get_Sector() != c2.get_Sector()
                                || c1.get_Sector() != c3.get_Sector()) {
                            continue;
                        }

                        double[] X = new double[3];
                        double[] Y = new double[3];
                        double[] Z = new double[3];
                        double[] errX = new double[3];
                        double[] errY = new double[3];

                        Z[0] = c1.get_Point().z();
                        Y[0] = c1.get_Point().y();
                        X[0] = c1.get_Point().x();
                        errX[0] = c1.get_PointErr().x();
                        errY[0] = c1.get_PointErr().y();
                        Z[1] = c2.get_Point().z();
                        Y[1] = c2.get_Point().y();
                        X[1] = c2.get_Point().x();
                        errX[1] = c2.get_PointErr().x();
                        errY[1] = c2.get_PointErr().y();
                        Z[2] = c3.get_Point().z();
                        Y[2] = c3.get_Point().y();
                        X[2] = c3.get_Point().x();
                        errX[2] = c3.get_PointErr().x();
                        errY[2] = c3.get_PointErr().y();

                        // Ignore point errors and assume the track vertex is close to the origin
                        TrajectoryParametriz qf1 = new TrajectoryParametriz();
                        qf1.evaluate(Z, X, errX, Y, errY);

                        Vector3D traj1 = new Vector3D(qf1.fitResult[3][0],
                                                      qf1.fitResult[4][0],
                                                      qf1.fitResult[5][0]);
                        Vector3D traj2 = new Vector3D(qf1.fitResult[3][1],
                                                      qf1.fitResult[4][1],
                                                      qf1.fitResult[5][1]);
                        Vector3D traj3 = new Vector3D(qf1.fitResult[3][2],
                                                      qf1.fitResult[4][2],
                                                      qf1.fitResult[5][2]);

                        double cosTh1 = traj1.dot(c1.get_Dir().toVector3D());
                        double cosTh2 = traj2.dot(c2.get_Dir().toVector3D());
                        double cosTh3 = traj3.dot(c3.get_Dir().toVector3D());

                        // Require that the cross direction estimate be in the  direction of the
                        //     trajectory
                        if (cosTh1 < Constants.TRACKDIRTOCROSSDIRCOSANGLE
                                || cosTh2 < Constants.TRACKDIRTOCROSSDIRCOSANGLE
                                || cosTh3 < Constants.TRACKDIRTOCROSSDIRCOSANGLE) {
                            continue;
                        }

                        double fitchsq = 0;

                        if (!c1.isPseudoCross) {
                            fitchsq += ((qf1.fitResult[1][0] - c1.get_Point().y()) /
                                        c1.get_PointErr().y()) *
                                       ((qf1.fitResult[1][0] - c1.get_Point().y()) /
                                        c1.get_PointErr().y());
                        }
                        if (!c2.isPseudoCross) {
                            fitchsq += ((qf1.fitResult[1][1] - c2.get_Point().y()) /
                                        c2.get_PointErr().y()) *
                                       ((qf1.fitResult[1][1] - c2.get_Point().y()) /
                                        c2.get_PointErr().y());
                        }
                        if (!c3.isPseudoCross) {
                            fitchsq += ((qf1.fitResult[1][2] - c3.get_Point().y()) /
                                        c3.get_PointErr().y()) *
                                       ((qf1.fitResult[1][2] - c3.get_Point().y()) /
                                        c3.get_PointErr().y());
                        }

                        // fit the  projection with a line
                        // The track is ~ constant in phi
                        LineFitter linefit = new LineFitter();
                        if (!linefit.fitStatus(X, Y, errX, errY, Z.length)) continue;

                        this.updateBFittedHits(c1, tab, DcDetector, tde, swimmer);
                        this.updateBFittedHits(c2, tab, DcDetector, tde, swimmer);
                        this.updateBFittedHits(c3, tab, DcDetector, tde, swimmer);

                        BaseCand bCand = new BaseCand();
                        bCand.CrossesOnTrack.clear();
                        bCand.CrossesOnTrack.add(c1);
                        bCand.CrossesOnTrack.add(c2);
                        bCand.CrossesOnTrack.add(c3);
                        bCand.Chisq = fitchsq;

                        if (bCand.Chisq < Constants.CROSSLISTSELECTQFMINCHSQ) trkCnds.add(bCand);
                    }
                }
            }
        }

        CrossList crossList = new CrossList();
        for (int i = 0; i < trkCnds.size(); i++) {
            crossList.add(i,trkCnds.get(i).CrossesOnTrack);
        }

        return crossList;
    }

    /**
     * Recalculates the direction of a cross
     * @param c1    the cross
     * @param slope
     */
    @SuppressWarnings("unused")
    private void RecalculateCrossDir(Cross c1, double slope) {
        double val_sl2 = c1.get_Segment2().get_fittedCluster().get_clusterLineFitSlope();
        double tanThX = val_sl2;
        double tanThY = slope;
        double uz = 1. / Math.sqrt(1 + tanThX*tanThX + tanThY*tanThY);
        double ux = uz*tanThX;
        double uy = uz*tanThY;

        c1.set_Dir(new Point3D(ux, uy, uz));
    }

    @SuppressWarnings("unused")
    private void RecalculateCrossDirErr(Cross c1,
                                        double slope,
                                        double slopeErr) {
        // Error calculation
        double val_sl2 = c1.get_Segment2().get_fittedCluster()
                                          .get_clusterLineFitSlope();
        double tanThX = val_sl2;
        double tanThY = slope;
        double del_tanThX = c1.get_Segment2().get_fittedCluster()
                                             .get_clusterLineFitSlopeErr();
        double del_tanThY = slopeErr;
        double uz = 1./Math.sqrt( 1 + tanThX*tanThX + tanThY*tanThY );
        double del_uz = uz*uz*uz*Math.sqrt(tanThX*tanThX*del_tanThX*del_tanThX +
                                           tanThY*tanThY*del_tanThY*del_tanThY);
        double del_ux = Math.sqrt(tanThX*tanThX*del_uz*del_uz +
                                  uz*uz*del_tanThX*del_tanThX);
        double del_uy = Math.sqrt(tanThY*tanThY*del_uz*del_uz +
                                  uz*uz*del_tanThY*del_tanThY);

        double err_xDir = del_ux;
        double err_yDir = del_uy;
        double err_zDir = del_uz;

        Point3D estimDirErr = new Point3D(err_xDir, err_yDir, err_zDir);

        c1.set_DirErr(estimDirErr);
    }

    private void recalcParsSegment(Segment _Segment, IndexedTable tab, DCGeant4Factory DcDetector,
            TimeToDistanceEstimator tde) {
        // Refit
        double cosTrkAngle = 1. / Math.sqrt(1. +
                _Segment.get_fittedCluster().get_clusterLineFitSlope() *
                _Segment.get_fittedCluster().get_clusterLineFitSlope());

        // Update the hits
        for (FittedHit fhit : _Segment.get_fittedCluster()) {
            fhit.updateHitPositionWithTime(cosTrkAngle, fhit.getB(), tab, DcDetector, tde);
        }

        cf.setFitArray(_Segment.get_fittedCluster(), true);
        cf.fit(_Segment.get_fittedCluster(), true);
        cosTrkAngle = 1. / Math.sqrt(1. + _Segment.get_fittedCluster().get_clusterLineFitSlope() *
                                          _Segment.get_fittedCluster().get_clusterLineFitSlope());

        for (FittedHit fhit : _Segment.get_fittedCluster()) {
            fhit.updateHitPositionWithTime(cosTrkAngle, fhit.getB(), tab, DcDetector, tde);
        }

        cf.setFitArray(_Segment.get_fittedCluster(), true);
        cf.fit(_Segment.get_fittedCluster(), true);

        cf.setResidualDerivedParams(_Segment.get_fittedCluster(), true, false, DcDetector);

        cf.setFitArray(_Segment.get_fittedCluster(), true);
        cf.fit(_Segment.get_fittedCluster(), false);

        cf.setSegmentLineParameters(_Segment.get_fittedCluster().get(0).get_Z(),
                                    _Segment.get_fittedCluster());
    }

    public List<List<Cross>> get_CrossesInSectors(List<Cross> crosses) {

        List<List<Cross>> CrossesBySectors = new ArrayList<List<Cross>>();
        for (int s = 0; s < 6; s++) CrossesBySectors.add(s, new ArrayList<Cross>());
        for (Cross cross : crosses) CrossesBySectors.get(cross.get_Sector()-1).add(cross);

        return CrossesBySectors;
    }

    /**
     * Updates B for the fitted hits in a cross.
     * @param c          cross
     * @param tab        table of constants
     * @param DcDetector detector geometry
     * @param tde        time-to-distance utility
     * @param swimmer
     * Updates the B-field information of the hits in the cross segments.
     */
    private void updateBFittedHits(Cross c, IndexedTable tab, DCGeant4Factory DcDetector,
            TimeToDistanceEstimator tde, Swim swimmer) {

        for (int i = 0; i < c.get_Segment1().size(); i++) {
            Point3D ref = c.get_Segment1().get(i).getCrossDirIntersWire();
            float[] result = new float[3];
            swimmer.Bfield(c.get_Sector(), ref.x(), ref.y(), ref.z(), result);
            c.get_Segment1().get(i).setB(Math.sqrt(result[0]*result[0] +
                                                   result[1]*result[1] +
                                                   result[2]*result[2]));
        }

        for (int i = 0; i < c.get_Segment2().size(); i++) {
            Point3D ref = c.get_Segment2().get(i).getCrossDirIntersWire();
            float[] result = new float[3];
            swimmer.Bfield(c.get_Sector(), ref.x(), ref.y(), ref.z(), result);
            c.get_Segment2().get(i).setB(Math.sqrt(result[0]*result[0] +
                                                   result[1]*result[1] +
                                                   result[2]*result[2]));
        }

        if (tde != null) {
            this.recalcParsSegment(c.get_Segment1(), tab, DcDetector, tde);
            this.recalcParsSegment(c.get_Segment2(), tab, DcDetector, tde);
        }

        // Remake cross
        c.set_CrossParams(DcDetector);
    }

    private class BaseCand {
        public double Chisq;
        public List<Cross> CrossesOnTrack = new ArrayList<Cross>();
    }

    private class TrajectoryParametriz {
        private double[][] fitResult = {{0.,0.,0.},
                                        {0.,0.,0.},
                                        {0.,0.,0.},
                                        {0.,0.,0.},
                                        {0.,0.,0.},
                                        {0.,0.,0.}};

        public double[] evaluate(double[] x, double[] y, double[] err, double[] y2, double[] err2) {

            LineFitter linefit = new LineFitter();
            linefit.fitStatus(x, y2, err, err, x.length);

            double[] ret = {0.,0.,0.};
            Matrix A = new Matrix(3,3);
            Matrix V = new Matrix(3,1);

            double sum1 = 0.0;
            double sum2 = 0.0;
            double sum3 = 0.0;
            double sum4 = 0.0;
            double sum5 = 0.0;
            double sum6 = 0.0;
            double sum7 = 0.0;
            double sum8 = 0.0;
            for (int i = 0; i < x.length; ++i) {
                double y1 = y[i];
                double x1 = x[i];
                double x2 = x1 * x1;
                double x3 = x2 * x1;
                double x4 = x2 * x2;
                double e2 = err[i] * err[i];

                sum1 += x4/e2;
                sum2 += x3/e2;
                sum3 += x2/e2;
                sum4 += x1/e2;
                sum5 += 1.0/e2;
                sum6 += y1 * x2/e2;
                sum7 += y1 * x1/e2;
                sum8 += y1/e2;
            }

            A.set(0, 0, sum1);
            A.set(0, 1, sum2);
            A.set(0, 2, sum3);
            A.set(1, 0, sum2);
            A.set(1, 1, sum3);
            A.set(1, 2, sum4);
            A.set(2, 0, sum3);
            A.set(2, 1, sum4);
            A.set(2, 2, sum5);
            V.set(0, 0, sum6);
            V.set(1, 0, sum7);
            V.set(2, 0, sum8);

            Matrix Ainv = A.inverse();
            Matrix X;

            try {
                X = Ainv.times(V);
                for (int i = 0; i < 3; ++i) {
                    ret[i] = X.get(i, 0);
                }
                for (int i = 0; i < x.length; i++) {

                    double tiltSysXterm = ret[0]*x[i]*x[i] + ret[1]*x[i] + ret[2];
                    double tiltSysYterm = linefit.getFit().slope()*x[i] +
                                          linefit.getFit().intercept();
                    double tiltSysZterm = x[i];

                    double dl  = 0.01;
                    double dQ  = 2.*ret[0]*x[i]*dl + ret[0] + ret[1]*dl;
                    double dL  = linefit.getFit().slope()*dl;
                    double Len = Math.sqrt(dl*dl + dQ*dQ + dL*dL) ;

                    double tiltSysdirXterm = dQ/Len;
                    double tiltSysdirYterm = dL/Len;
                    double tiltSysdirZterm = dl/Len;

                    fitResult[0][i] = tiltSysXterm;
                    fitResult[1][i] = tiltSysYterm;
                    fitResult[2][i] = tiltSysZterm;
                    fitResult[3][i] = tiltSysdirXterm;
                    fitResult[4][i] = tiltSysdirYterm;
                    fitResult[5][i] = tiltSysdirZterm;
                }
            }
            catch (ArithmeticException e) {
                e.printStackTrace();
            }
            return(ret);
        }
    }
}
