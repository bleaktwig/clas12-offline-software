package org.jlab.rec.dc.cluster;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.jlab.geom.prim.Point3D;
import org.jlab.geom.prim.Vector3D;
import org.jlab.geom.prim.Line3D;
import org.jlab.detector.geant4.v2.DCGeant4Factory;

import org.jlab.rec.dc.track.fit.basefit.LineFitPars;
import org.jlab.rec.dc.track.fit.basefit.LineFitter;
import org.jlab.rec.dc.hit.FittedHit;

/**
 * Fits a cluster to a line.
 *
 * @author ziegler
 */
public class ClusterFitter {

    private LineFitPars FitPars;
    private List<ArrayList<Double>> FitArray = new ArrayList<ArrayList<Double>>();
    private List<Double> x = new ArrayList<Double>();
    private List<Double> y = new ArrayList<Double>();
    private List<Double> ex = new ArrayList<Double>();
    private List<Double> ey = new ArrayList<Double>();
    // Math.cos(Math.toRadians(6.)) = 0.9945218953682733
    private double stereo = 0.9945218953682733;

    private String CoordinateSystem; // LC = local, TSC = tilted Sector
                                     // LC:  local coordinate grid, Delta_z = 1
                                     // TSC: local tilted coordinate system Delta_z, ~ cell size

    public ClusterFitter() {}

    /**
     * NOTE: Lacks description
     * @param clus   NOTE: Lacks description
     * @param system NOTE: Lacks description
     */
    public void SetFitArray(FittedCluster clus, String system) {

        Collections.sort(clus);
        for (int i = 0; i < FitArray.size(); i++) FitArray.get(i).clear();

        for (int i = 0; i < clus.size(); i++) {
            if (system.equals("LC")) {
                CoordinateSystem = "LC";
                x.add(i, clus.get(i).get_lX());
                ex.add(i, (double) 0);
                y.add(i, clus.get(i).get_lY());
                ey.add(i, (double) 1);
            }
            else if (system.equals("TSC")) {
                CoordinateSystem = "TSC";
                x.add(i, clus.get(i).get_Z());
                ex.add(i, (double) 0);
                y.add(i, clus.get(i).get_X());
                ey.add(i, clus.get(i).get_DocaErr() / stereo);
            }
        }
        FitArray.add((ArrayList<Double>) x);
        FitArray.add((ArrayList<Double>) ex);
        FitArray.add((ArrayList<Double>) y);
        FitArray.add((ArrayList<Double>) ey);
    }

    /**
     * Fits the cluster to a line if it hasn't been done yet.
     * @param clus fitted cluster
     * @param SaveFitPars boolean to save the fit parameters
     */
    public void Fit (FittedCluster clus, boolean SaveFitPars) {
        if (FitArray != null) {
            LineFitter linefit = new LineFitter();
            boolean linefitstatusOK = linefit.fitStatus(FitArray.get(0), FitArray.get(2),
                                                        FitArray.get(1), FitArray.get(3),
                                                        FitArray.get(0).size());

            if (linefitstatusOK) FitPars = linefit.getFit();
            else FitPars = null;

            if (SaveFitPars) this.SetClusterFitParameters(clus);
        }
    }

    /**
     * Saves the fit parameters into a fitted cluster.
     * @param clus fitted cluster
     */
    public void SetClusterFitParameters(FittedCluster clus) {
        if (FitPars == null) return;

        clus.set_clusterLineFitSlope(FitPars.slope());
        clus.set_clusterLineFitSlopeErr(FitPars.slopeErr());
        clus.set_clusterLineFitIntercept(FitPars.intercept());
        clus.set_clusterLineFitInterceptErr(FitPars.interceptErr());
        clus.set_clusterLineFitSlIntCov(FitPars.SlopeIntercCov());
        clus.set_fitProb(FitPars.getProb());
        clus.set_Chisq(FitPars.chisq());
    }

    /**
     * NOTE: Method lacks description
     * @param x0   local x in the tilted sector coordinate system (in cm)
     * @param clus fitted cluster
     */
    public void SetSegmentLineParameters(double x0, FittedCluster clus) {
        if (FitPars == null) return;

        Point3D pointOnLine    = get_PointOnLine(x0, FitPars.slope(), FitPars.intercept());
        Point3D dirOfLine      = get_DirOnLine(FitPars.slope(), FitPars.intercept());
        Point3D pointOnLineErr = get_PointOnLine(x0, FitPars.slopeErr(), FitPars.interceptErr());
        Point3D dirOfLineErr   = get_DirOnLine(FitPars.slopeErr(), FitPars.interceptErr());

        clus.set_clusLine(new Line3D(pointOnLine, dirOfLine));
        clus.set_clusLineErr(new Line3D(pointOnLineErr, dirOfLineErr));
        clus.set_clusterLineFitSlopeMP(FitPars.slope());
        clus.set_clusterLineFitSlopeErrMP(FitPars.slopeErr());
        clus.set_clusterLineFitInterceptMP(FitPars.intercept());
        clus.set_clusterLineFitInterceptErrMP(FitPars.interceptErr());
    }

    /**
     * NOTE: Lacks description
     * @param clus             fitted cluster
     * @param calcTimeResidual boolean to compute the time residuals (in cm)
     * @param resetLRAmbig     boolean to reset the LR ambiguity
     * @param DcDetector       DC detector geometry
     */
    public void SetResidualDerivedParams(FittedCluster clus,
                                         boolean calcTimeResidual,
                                         boolean resetLRAmbig,
                                         DCGeant4Factory DcDetector) {

        if (FitPars == null || FitArray == null) return;

        // 6 = nb of layers, 3 number of criteria
        int[][] statusArray     = new int[3][6];
        int[] nHitsInLyr        = new int[6];
        int[] nLRresolvedInLyr  = new int[6];
        int[] nPassingResAccCut = new int[6];

        for (int i = 0; i < clus.size(); i++) {
            nHitsInLyr[clus.get(i).get_Layer() - 1]++;

            double residual = (FitArray.get(2).get(i) - FitPars.slope() * FitArray.get(0).get(i) -
                               FitPars.intercept());
            clus.get(i).set_Residual(residual);

            double xWire = DcDetector.getWireMidpoint(clus.get(i).get_Sector() - 1,
                                                      clus.get(i).get_Superlayer() - 1,
                                                      clus.get(i).get_Layer() - 1,
                                                      clus.get(i).get_Wire() - 1).x;
            double zWire = DcDetector.getWireMidpoint(clus.get(i).get_Sector() - 1,
                                                      clus.get(i).get_Superlayer() - 1,
                                                      clus.get(i).get_Layer() - 1,
                                                      clus.get(i).get_Wire() - 1).z;

            Line3D FitLine = new Line3D();
            Point3D pointOnTrk = new Point3D(FitArray.get(0).get(0),
                    FitPars.slope() * FitArray.get(0).get(0) + FitPars.intercept(),
                    0);

            Vector3D trkDir = new Vector3D(1, FitPars.slope(), 0);
            trkDir.unit();
            FitLine.set(pointOnTrk, trkDir);
            Point3D Wire = new Point3D(zWire, xWire, 0);

            double trkDocaMP = FitLine.distance(Wire).length();
            double trkDoca   = trkDocaMP * stereo;

            clus.get(i).set_ClusFitDoca(trkDoca);

            if (Math.abs(residual) < 0.350) {
                // less than the average resolution
                nPassingResAccCut[clus.get(i).get_Layer() - 1]++;
            }

            if (clus.get(i).get_LeftRightAmb() == 0) {
                if      (residual < 0) clus.get(i).set_LeftRightAmb(1);
                else if (residual > 0) clus.get(i).set_LeftRightAmb(-1);
            }

            if (resetLRAmbig) {
                if ((CoordinateSystem.equals("LC") && Math.abs(residual) < 0.01)
                        || (CoordinateSystem.equals("LTS")
                        && clus.get(i).get_Doca() / clus.get(i).get_CellSize() < 0.4)) {
                    // DOCA is required to be larger than 40% of cell size for hit-based tracking LR
                    //     assignment
                    clus.get(i).set_LeftRightAmb(0);
                }
            }

            if (clus.get(i).get_LeftRightAmb() != 0) {
                nLRresolvedInLyr[clus.get(i).get_Layer() - 1]++;
            }

            if (calcTimeResidual == true) {
                double timeResidual = Math.abs(FitPars.slope() * FitArray.get(0).get(i)
                                    + FitPars.intercept() - xWire)
                                    - Math.abs(FitArray.get(2).get(i) - xWire);
                clus.get(i).set_TimeResidual(timeResidual);
            }
        }

        for (int i = 0; i < 6; i++) {
            statusArray[0][i] = nHitsInLyr[i];
            if (nHitsInLyr[i] > 0) {
                statusArray[1][i] = nLRresolvedInLyr[i]  / nHitsInLyr[i];
                statusArray[2][i] = nPassingResAccCut[i] / nHitsInLyr[i];
            }
        }

        clus.set_Status(statusArray);
    }

    /**
     * Selects the fitted cluster with the best fit chi^2.
     * @param clusters fitted cluster
     * @param system   coordinate system in which the fit is performed
     * @return         the selected fitted cluster
     */
    public FittedCluster BestClusterSelector(List<FittedCluster> clusters, String system) {

        FittedCluster BestCluster = null;
        double bestChisq = 999999999.;

        for (FittedCluster clusCand : clusters) {
            if (isBrickWall(clusCand)) {
                int LRSum = 0;
                for (FittedHit hit : clusCand) LRSum += hit.get_LeftRightAmb();
                if (LRSum != 0) continue;
            }

            // Set the array of measurements according to the system used in the analysis
            SetFitArray(clusCand, system);
            // Do the fit and get the chisq
            Fit(clusCand, true);
            if (FitPars == null) continue;
            double chisq = FitPars.chisq();

            if (chisq < bestChisq) {
                bestChisq = chisq;
                BestCluster = clusCand;
            }
        }
        return BestCluster;
    }

    /**
     * @param the_slope  the cluster fitted-line slope
     * @param the_interc the fitted-line intercept
     * @return           the cluster fitted-line unit direction vector
     */
    private Point3D get_DirOnLine(double the_slope, double the_interc) {
        return new Point3D(the_slope / Math.sqrt(1. + the_slope * the_slope), 0,
                           1. / Math.sqrt(1. + the_slope * the_slope));
    }

    /**
     * @param d          the point z coordinate
     * @param the_slope  the cluster fitted-line slope
     * @param the_interc the fitted-line intercept
     * @return           point (the_slope*d+the_interc,0,d) on the fitted cluster line
     */
    private Point3D get_PointOnLine(double d, double the_slope, double the_interc) {
        return new Point3D(the_slope * d + the_interc, 0, d);
    }

    /**
     * NOTE: Method lacks description
     * @param clusCand fitted cluster
     * @return wire pattern in the cluster
     */
    private boolean isBrickWall(FittedCluster clusCand) {
        boolean isBW = true;
        int sumWireNum = 0;
        if (clusCand.size() != 6) isBW = false;

        for (FittedHit hit : clusCand) sumWireNum+=hit.get_Wire();
        for (FittedHit hit : clusCand) {
            if (hit.get_Wire() * clusCand.size() != sumWireNum) isBW = false;
        }
        return isBW;
    }
}
