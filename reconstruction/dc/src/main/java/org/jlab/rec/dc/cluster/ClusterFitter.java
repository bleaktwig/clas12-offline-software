package org.jlab.rec.dc.cluster;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.jlab.detector.geant4.v2.DCGeant4Factory;
import org.jlab.geom.prim.Point3D;
import org.jlab.geom.prim.Vector3D;
import org.jlab.geom.prim.Line3D;
import org.jlab.rec.dc.hit.FittedHit;
import org.jlab.rec.dc.track.fit.basefit.LineFitPars;
import org.jlab.rec.dc.track.fit.basefit.LineFitter;

/**
 * Fits a cluster to a line.
 *
 * @author ziegler
 */
public class ClusterFitter {

    private LineFitPars fitPars;
    private List<ArrayList<Double>> fitArray = new ArrayList<ArrayList<Double>>();
    private List<Double> x = new ArrayList<Double>();
    private List<Double> y = new ArrayList<Double>();
    private List<Double> ex = new ArrayList<Double>();
    private List<Double> ey = new ArrayList<Double>();
    // Math.cos(Math.toRadians(6.)) = 0.9945218953682733
    private double stereo = 0.9945218953682733;

    private boolean coordinateSystem; // LC = local, TSC = tilted Sector
                                      // false - LC:  local coordinate grid, Delta_z = 1
                                      // true  - TSC: local tilted coordinate system Delta_z, ~ cell size

    public ClusterFitter() {}

    /**
     * @param system decides which system is being used.
     *               false: LC, local coordinate grid
     *               true:  TSC, local tilted coordinate system
     */
    public void setFitArray(FittedCluster clus, boolean system) {

        Collections.sort(clus);
        for (int i = 0; i < fitArray.size(); i++) fitArray.get(i).clear();

        this.coordinateSystem = system;
        if (!system) {
            for (int i = 0; i < clus.size(); i++) {
                x.add(i, clus.get(i).get_lX());
                ex.add(i, (double) 0);
                y.add(i, clus.get(i).get_lY());
                ey.add(i, (double) 1);
            }
        }
        else {
            for (int i = 0; i < clus.size(); i++) {
                x.add(i, clus.get(i).get_Z());
                ex.add(i, (double) 0);
                y.add(i, clus.get(i).get_X());
                ey.add(i, clus.get(i).get_DocaErr() / stereo);
            }
        }

        fitArray.add((ArrayList<Double>) x);
        fitArray.add((ArrayList<Double>) ex);
        fitArray.add((ArrayList<Double>) y);
        fitArray.add((ArrayList<Double>) ey);
    }

    /**
     * Fits the cluster to a line if it hasn't been done yet.
     * @param clus        fitted cluster
     * @param saveFitPars boolean to save the fit parameters
     */
    public void fit (FittedCluster clus, boolean saveFitPars) {

        if (fitArray == null) return;
        LineFitter lineFit = new LineFitter();
        boolean lineFitStatusOK = lineFit.fitStatus(fitArray.get(0), fitArray.get(2),
                                                    fitArray.get(1), fitArray.get(3),
                                                    fitArray.get(0).size());

        if (lineFitStatusOK) fitPars = lineFit.getFit();
        else                 fitPars = null;

        if (saveFitPars) this.setClusterFitParameters(clus);
    }

    /**
     * Saves the fit parameters into a fitted cluster.
     * @param clus fitted cluster
     */
    private void setClusterFitParameters(FittedCluster clus) {
        if (fitPars == null) return;

        clus.set_clusterLineFitSlope(fitPars.slope());
        clus.set_clusterLineFitSlopeErr(fitPars.slopeErr());
        clus.set_clusterLineFitIntercept(fitPars.intercept());
        clus.set_clusterLineFitInterceptErr(fitPars.interceptErr());
        clus.set_clusterLineFitSlIntCov(fitPars.SlopeIntercCov());
        clus.set_fitProb(fitPars.getProb());
        clus.set_Chisq(fitPars.chisq());
    }

    /**
     * @param x0   local x in the tilted sector coordinate system (in cm)
     * @param clus fitted cluster
     */
    public void setSegmentLineParameters(double x0, FittedCluster clus) {
        if (fitPars == null) return;

        Point3D pointOnLine    = get_PointOnLine(x0, fitPars.slope(), fitPars.intercept());
        Point3D dirOfLine      = get_DirOnLine(fitPars.slope(), fitPars.intercept());
        Point3D pointOnLineErr = get_PointOnLine(x0, fitPars.slopeErr(), fitPars.interceptErr());
        Point3D dirOfLineErr   = get_DirOnLine(fitPars.slopeErr(), fitPars.interceptErr());

        clus.set_clusLine(new Line3D(pointOnLine, dirOfLine));
        clus.set_clusLineErr(new Line3D(pointOnLineErr, dirOfLineErr));
        clus.set_clusterLineFitSlopeMP(fitPars.slope());
        clus.set_clusterLineFitSlopeErrMP(fitPars.slopeErr());
        clus.set_clusterLineFitInterceptMP(fitPars.intercept());
        clus.set_clusterLineFitInterceptErrMP(fitPars.interceptErr());
    }

    /**
     * @param clus             fitted cluster
     * @param calcTimeResidual boolean to compute the time residuals (in cm)
     * @param resetLRAmbig     boolean to reset the LR ambiguity
     * @param DcDetector       DC detector geometry
     */
    public void setResidualDerivedParams(FittedCluster clus, boolean calcTimeResidual,
            boolean resetLRAmbig, DCGeant4Factory DcDetector) {

        if (fitPars == null || fitArray == null) return;

        // 6 = nb of layers, 3 number of criteria
        int[][] statusArray     = new int[3][6];
        int[] nHitsInLyr        = new int[6];
        int[] nLRresolvedInLyr  = new int[6];
        int[] nPassingResAccCut = new int[6];

        for (int i = 0; i < clus.size(); i++) {
            nHitsInLyr[clus.get(i).get_Layer() - 1]++;

            clus.get(i).set_Residual(fitArray.get(2).get(i) - fitPars.slope()
                                     * fitArray.get(0).get(i) - fitPars.intercept());

            double xWire = DcDetector.getWireMidpoint(clus.get(i).get_Sector() - 1,
                                                      clus.get(i).get_Superlayer() - 1,
                                                      clus.get(i).get_Layer() - 1,
                                                      clus.get(i).get_Wire() - 1).x;
            double zWire = DcDetector.getWireMidpoint(clus.get(i).get_Sector() - 1,
                                                      clus.get(i).get_Superlayer() - 1,
                                                      clus.get(i).get_Layer() - 1,
                                                      clus.get(i).get_Wire() - 1).z;

            Point3D pointOnTrk = new Point3D(fitArray.get(0).get(0),
                    fitPars.slope() * fitArray.get(0).get(0) + fitPars.intercept(), 0);

            Vector3D trkDir = new Vector3D(1, fitPars.slope(), 0);
            trkDir.unit();

            Line3D fitLine = new Line3D();
            fitLine.set(pointOnTrk, trkDir);

            Point3D wire = new Point3D(zWire, xWire, 0);
            double trkDocaMP = fitLine.distance(wire).length();
            clus.get(i).set_ClusFitDoca(trkDocaMP * stereo);

            if (Math.abs(clus.get(i).get_Residual()) < 0.350) {
                // less than the average resolution
                nPassingResAccCut[clus.get(i).get_Layer() - 1]++;
            }

            if (clus.get(i).get_LeftRightAmb() == 0) {
                if      (clus.get(i).get_Residual() < 0) clus.get(i).set_LeftRightAmb(1);
                else if (clus.get(i).get_Residual() > 0) clus.get(i).set_LeftRightAmb(-1);
            }

            if (resetLRAmbig && ((!coordinateSystem && Math.abs(clus.get(i).get_Residual())<0.01)
                    || (coordinateSystem && clus.get(i).get_Doca()/clus.get(i).get_CellSize()<0.4)))
                // DOCA is required to be larger than 40% of cell size for hit-based tracking LR
                //     assignment
                clus.get(i).set_LeftRightAmb(0);

            if (clus.get(i).get_LeftRightAmb() != 0)
                nLRresolvedInLyr[clus.get(i).get_Layer() - 1]++;

            if (calcTimeResidual) {
                double timeResidual = Math.abs(fitPars.slope() * fitArray.get(0).get(i)
                                    + fitPars.intercept() - xWire)
                                    - Math.abs(fitArray.get(2).get(i) - xWire);
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
    public FittedCluster BestClusterSelector(List<FittedCluster> clusters, boolean system) {

        FittedCluster BestCluster = null;
        double bestChisq = 999999999.;

        for (FittedCluster clusCand : clusters) {
            if (isBrickWall(clusCand)) {
                int LRSum = 0;
                for (FittedHit hit : clusCand) LRSum += hit.get_LeftRightAmb();
                if (LRSum != 0) continue;
            }

            // Set the array of measurements according to the system used in the analysis
            setFitArray(clusCand, system);
            // Do the fit and get the chisq
            fit(clusCand, true);
            if (fitPars == null) continue;
            double chisq = fitPars.chisq();

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
