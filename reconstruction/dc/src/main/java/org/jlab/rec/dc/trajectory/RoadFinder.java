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
    public RoadFinder() {}

    /**
     * Finds a list of pseudo-segments
     * @param segs       list of segments
     * @param DcDetector DC detector utility
     * @return           list of segments corresponding to pseudo-segments
     */
    public List<Road> findRoads(List<Segment> segs, DCGeant4Factory DcDetector) {

        // Initialize the lists
        List<Road> roads = new ArrayList<Road>();

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

                            if (sLyr.size() < 3 || !this.fitRoad(sLyr, DcDetector) || qf.chi2 >= 150
                                    || qf.chi2 == 0) {
                                continue;
                            }
                            // Road is good --> pass w.out looking for missing segment
                            sLyr.id = roadId;
                            sLyr.a  = qf.a;
                            roads.add(sLyr);
                            roadId++;
                        }
                    }
                }
            }
        }

        return roads;
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

        cf.setFitArray(fPseudoCluster, true);
        cf.fit(fPseudoCluster, true);

        cf.setSegmentLineParameters(fPseudoCluster.get(0).get_Z(), fPseudoCluster);
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
        cf.setFitArray(fPseudoCluster, true);
        cf.fit(fPseudoCluster, true);

        cf.setSegmentLineParameters(fPseudoCluster.get(0).get_Z(), fPseudoCluster);
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
            double[] sum = {0, 0, 0, 0, 0, 0, 0, 0};

            for (int i = 0; i < x.length; ++i) {
                double y1 = y[i];
                double x1 = x[i];
                double x2 = x1 * x1;
                double x3 = x2 * x1;
                double x4 = x2 * x2;
                double e2 = err[i] * err[i];

                sum[0] += x4/e2;
                sum[1] += x3/e2;
                sum[2] += x2/e2;
                sum[3] += x1/e2;
                sum[4] += 1.0/e2;
                sum[5] += y1 * x2/e2;
                sum[6] += y1 * x1/e2;
                sum[7] += y1/e2;
            }

            // NOTE: The matrix inversion is hardcoded so as to accelerate the method
            double g1 = sum[2]*sum[2];
            double g2 = sum[2]*sum[3];
            double g3 = sum[1]*sum[4];
            double g4 = sum[0]*sum[3];
            double g5 = sum[2]*sum[4];

            double i1 = g2 - g3;
            double i2 = sum[1]*sum[3] - g1;
            double i3 = sum[1]*sum[2] - g4;

            double[] ret = {(g5 - sum[3]*sum[3]) * sum[5] + i1 * sum[6] + i2 * sum[7],
                            i1 * sum[5] + (sum[0]*sum[4] - g1) * sum[6] + i3 * sum[7],
                            i2 * sum[5] + i3 * sum[6] + (sum[0]*sum[2] - sum[1]*sum[1]) * sum[7]};

            double div = -g1*sum[2] + 2*sum[1]*g2 + sum[0]*g5 - sum[3]*g4 - sum[1]*g3;

            for (int i = 0; i < ret.length; ++i) ret[i] /= div;

            double _chi2 = 0;
            for (int i = 0; i < x.length; ++i) {
                double tiltSysXterm = ret[0]*x[i]*x[i] + ret[1]*x[i] + ret[2];
                _chi2 += (tiltSysXterm-y[i]) * (tiltSysXterm-y[i]) / (err[i]*err[i]);
            }
            this.chi2 = _chi2;
            this.NDF = x.length - 3;

            a = ret;
        }
    }
}
