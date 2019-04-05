package org.jlab.rec.dc.trajectory;

import java.util.ArrayList;
import java.util.List;
import Jama.Matrix;

import org.jlab.detector.geant4.v2.DCGeant4Factory;

import org.jlab.rec.dc.hit.FittedHit;
import org.jlab.rec.dc.cluster.Cluster;
import org.jlab.rec.dc.cluster.FittedCluster;
import org.jlab.rec.dc.cluster.ClusterFitter;
import org.jlab.rec.dc.segment.Segment;

/**
 * @author ziegler
 */
public class RoadFinder  {

    private SegmentTrajectory segTrj = new SegmentTrajectory();

    private ClusterFitter cf = new ClusterFitter();
    QuadraticFit qf = new QuadraticFit();
    public RoadFinder() {
    }

    /**
     * Finds a list of pseudo-segments
     * @param segs       list of segments
     * @param DcDetector DC detector utility
     * @return           list of segments corresponding to pseudo-segments
     */
    public List<Road> findRoads(List<Segment> segs, DCGeant4Factory DcDetector) {

        // Initialize the lists
        List<Road> Roads = new ArrayList<Road>();

        List<ArrayList<ArrayList<Segment>>> superLayerLists =
                new ArrayList<ArrayList<ArrayList<Segment>>>();
        for (int sec = 0; sec < 6; sec++) {
            ArrayList<ArrayList<Segment>> sLyrs = new ArrayList<ArrayList<Segment>>();
            ArrayList<ArrayList<ArrayList<Segment>>> rLyrs =
                    new ArrayList<ArrayList<ArrayList<Segment>>>();

            for (int sly = 0; sly < 6; sly++) sLyrs.add(new ArrayList<Segment>());
            superLayerLists.add(sLyrs);
        }

        // Make an array sorted by sector, superlayers
        for (Segment seg : segs) {
            if (seg.isOnTrack == false) {
                superLayerLists.get(seg.get_Sector() - 1)
                               .get(seg.get_Superlayer()-1)
                               .add((Segment) seg.clone());
            }
        }

        for (int sec = 0; sec < 6; sec++)  {
            for (int sly = 0; sly < 6; sly++) {
                Segment blank = new Segment(new FittedCluster(new Cluster(sec + 1, sly + 1, -1)));
                blank.set_Id(-10);
                superLayerLists.get(sec).get(sly).add(blank);
            }
        }

        int roadId = 1;
        for (int sec = 0; sec < 6; sec++) {
            for (int j = 0; j < 2; j++) {
                for (Segment s1 : superLayerLists.get(sec).get(0 + j)) {
                    for (Segment s2 : superLayerLists.get(sec).get(2 + j)) {
                        for (Segment s3 : superLayerLists.get(sec).get(4 + j)) {
                            Road sLyr = new Road();
                            if (s1.get_Id() != -10) sLyr.add(s1);
                            if (s2.get_Id() != -10) sLyr.add(s2);
                            if (s3.get_Id() != -10) sLyr.add(s3);

                            if (sLyr.size() < 3
                                    || !this.fitRoad(sLyr, DcDetector)
                                    || qf.chi2 >= 150
                                    || qf.chi2 == 0) {
                                continue;
                            }
                            // Road is good --> pass w.out looking for missing segment
                            sLyr.id = roadId;
                            sLyr.a  = qf.a;
                            Roads.add(sLyr);
                            roadId++;
                        }
                    }
                }
            }
        }

        return Roads;
    }

    public Segment findRoadMissingSegment(List<Segment> segList,
                                          DCGeant4Factory DcDetector,
                                          double[] a)  {

        if (segList.size() >= 3) return null;
        // else
        Segment pseudoSeg = null;
        // Make pseudo-segment for missing segment
        // Find missing segment superlayer
        int s1 = (segList.get(0).get_Superlayer() - (segList.get(0).get_Superlayer()+1)%2-1)/2;
        int s2 = (segList.get(1).get_Superlayer() - (segList.get(1).get_Superlayer()+1)%2-1)/2;
        int smiss = -1;
        if (s1 == 0) {
            if (s2 == 1) smiss = 2;
            if (s2 == 2) smiss = 1;
        } else {
            smiss = 0;
        }

        // The missing superlayer
        int slyr = (segList.get(0).get_Superlayer() + 1)%2 + 2*smiss + 1;
        if (slyr < 1 || slyr > 6) return null;

        // Make the missing segment
        Cluster pseudoCluster = new Cluster(segList.get(0).get_Sector(), slyr, -1);
        FittedCluster fPseudoCluster = new FittedCluster(pseudoCluster);
        for (int l = 0; l < 6; l++) {
            int layer = l + 1;
            double z = DcDetector.getWireMidpoint(segList.get(0).get_Sector() - 1,
                                                  slyr - 1, layer - 1, 0).z;
            double trkX = a[0]*z*z + a[1]*z + a[2];
            int calcWire = segTrj.getWireOnTrajectory(segList.get(0).get_Sector(),
                                                      slyr, layer, trkX, DcDetector);
            FittedHit pseudoHit = new FittedHit(segList.get(0).get_Sector(),
                                                slyr, layer, calcWire, 0, -1);
            // Estimate the error on the hit as the cellSize/sqrt(12)
            pseudoHit.calc_CellSize(DcDetector);
            // Math.sqrt(12.)) = 3.4641016151377544
            pseudoHit.set_DocaErr(pseudoHit.get_CellSize() / 3.4641016151377544);
            // Update the hit position estimate and add to the pseudo-cluster
            pseudoHit.updateHitPosition(DcDetector);
            fPseudoCluster.add(pseudoHit);
        }

        cf.SetFitArray(fPseudoCluster, "TSC");
        cf.Fit(fPseudoCluster, true);

        cf.SetSegmentLineParameters(fPseudoCluster.get(0).get_Z(), fPseudoCluster);
        pseudoSeg = new Segment(fPseudoCluster);
        pseudoSeg.set_fitPlane(DcDetector);
        return pseudoSeg;
    }

    /**
     * Redo a segment's fit (not used).
     * @param pseudoSeg  segment to be re-fitted
     * @param segList    list of segments
     * @param DcDetector DC detector geometry
     * @return           the new segment
     */
    private Segment reFit(Segment pseudoSeg,
                          ArrayList<Segment> segList,
                          DCGeant4Factory DcDetector ) {

        qf.init();
        this.fitRoad(segList, DcDetector);

        Cluster pseudoCluster = new Cluster(segList.get(0).get_Sector(),
                                            pseudoSeg.get_Superlayer(), -1);
        FittedCluster fPseudoCluster = new FittedCluster(pseudoCluster);

        for (int l = 0; l < 6; l++) {
            int layer = l + 1;
            double z = DcDetector.getWireMidpoint(pseudoSeg.get_Sector() - 1,
                                                  pseudoSeg.get_Superlayer() - 1,
                                                  layer - 1, 0).z;

            double trkX = qf.a[0]*z*z + qf.a[1]*z + qf.a[2];
            // Math.cos(Math.toRadians(6.)) = 0.9945218953682733
            double delta = (trkX-pseudoSeg.get(l).get_X())
                         / pseudoSeg.get(l).get_CellSize() / 0.9945218953682733;
            int calcWire = segTrj.getWireOnTrajectory(pseudoSeg.get_Sector(),
                                                      pseudoSeg.get_Superlayer(),
                                                      layer, trkX, DcDetector);

            FittedHit pseudoHit = new FittedHit(segList.get(0).get_Sector(),
                                                pseudoSeg.get_Superlayer(),
                                                layer, calcWire, 0, -1);

            // Math.sqrt(12.) / Math.cos(Math.toRadians(6.)) = 3.483182855270362
            pseudoHit.set_DocaErr(pseudoHit.get_CellSize() / 3.483182855270362);
            pseudoHit.updateHitPosition(DcDetector);
            fPseudoCluster.add(pseudoHit);
        }
        cf.SetFitArray(fPseudoCluster, "TSC");
        cf.Fit(fPseudoCluster, true);

        cf.SetSegmentLineParameters(fPseudoCluster.get(0).get_Z(), fPseudoCluster) ;
        Segment pseudoSeg1 = new Segment(fPseudoCluster);

        pseudoSeg1.set_fitPlane(DcDetector);

        return pseudoSeg1;
    }

    private boolean fitRoad(ArrayList<Segment> segList, DCGeant4Factory DcDetector) {

        qf.init();
        int NbHits = 0;

        if (segList.size() < 2) return false;
        for (Segment s : segList) NbHits += s.size();

        double[] X    = new double[NbHits];
        double[] Z    = new double[NbHits];
        double[] errX = new double[NbHits];

        int hitno = 0;

        // Math.sqrt(12.) / Math.cos(Math.toRadians(6.)) = 3.483182855270362
        for (Segment s : segList) {
            for (int j = 0; j < s.size(); j++) {
                X[hitno]    = s.get(j).get_X();
                Z[hitno]    = s.get(j).get_Z();
                errX[hitno] = s.get(j).get_CellSize() / 3.483182855270362;
                hitno++;
            }
        }

        qf.evaluate(Z, X, errX);

        double WChi2 = 0;
        for (Segment s : segList) {
            for (FittedHit h : s) {
                double trkX = qf.a[0]*h.get_Z()*h.get_Z() + qf.a[1]*h.get_Z() + qf.a[2];
                int calcWire = segTrj.getWireOnTrajectory(h.get_Sector(), h.get_Superlayer(),
                                                          h.get_Layer(), trkX, DcDetector);
                WChi2 += (h.get_Wire() - calcWire) * (h.get_Wire() - calcWire);
            }
        }

        // Pass if normalized chi2 is less than 1
        return WChi2 / qf.NDF <= 1;
    }

    /**
     * Finds the superlayer in a region in which there should be a match to make a cross; i.e. for
     * a segment in superlayer 1 there should be a matched segment in superlayer 2.
     * @param superlayer segment superlayer
     * @return           the superlayer
     */
    private int SuperlayerInWhichToSearchMatchingSeg(int superlayer) {
        if (superlayer % 2 == 0) return superlayer - 1; // Even layer
        else                     return superlayer + 1; // Odd layer
    }

    /** Quadratic fitting class to handle trajectory approximation. */
    private class QuadraticFit {
        public double chi2;
        public double NDF;
        public double[] a;

        public void init() {
            chi2 = 0;
            NDF  = 0;
            a    = new double[3];
        }

        public void evaluate(double[] x, double[] y, double[] err) {

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

            // TODO: Hardcoding the inverse so it doesn't need to be calculated
            // double aux1 = sum3*sum4 + sum2*sum5;
            // double aux2 = sum3*sum3;
            // double aux3 = sum2*sum4 - aux2;
            // double aux4 = sum2*sum3 - sum1*sum4;
            //
            // double[] ret = {(sum3*sum5 - sum4*sum4)*sum6 + aux1*sum7 + aux3*sum8,
            //                 aux1*sum6 + (sum1*sum5 - aux2)*sum7 + aux4*sum8,
            //                 aux3*sum6 + aux4*sum7 + (sum1*sum3 - sum2*sum2)*sum8};
            //
            // double _chi2 = 0;
            // for (int i = 0; i < x.length; ++i) {
            //     double tiltSysXterm = ret[0]*x[i]*x[i] + ret[1]*x[i] + ret[2];
            //     _chi2 += (tiltSysXterm-y[i]) * (tiltSysXterm-y[i]) / (err[i]*err[i]);
            // }
            // this.chi2 = _chi2;
            // this.NDF = x.length - 3;
            //
            // a = ret;

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
                for (int i = 0; i < 3; ++i) ret[i] = X.get(i, 0);

                double _chi2 = 0;
                for (int i = 0; i < x.length; i++) {
                    double tiltSysXterm = ret[0]*x[i]*x[i] + ret[1]*x[i] + ret[2];
                    _chi2 += (tiltSysXterm-y[i]) * (tiltSysXterm-y[i]) / (err[i]*err[i]);
                }
                this.chi2 = _chi2;
                this.NDF  = x.length -3;
            } catch (ArithmeticException e) {
            }
            a = ret;
        }
    }
}
