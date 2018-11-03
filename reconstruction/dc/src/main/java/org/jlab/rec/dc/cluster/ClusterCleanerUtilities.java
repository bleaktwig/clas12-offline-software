package org.jlab.rec.dc.cluster;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.jlab.utils.groups.IndexedTable;
import org.jlab.detector.geant4.v2.DCGeant4Factory;

import org.jlab.rec.dc.Constants;
import org.jlab.rec.dc.timetodistance.TimeToDistanceEstimator;
import org.jlab.rec.dc.hit.Hit;
import org.jlab.rec.dc.hit.FittedHit;

// NOTE: Maybe changing the Math methods here to FastMath would be practical?
// NOTE: Lacks JavaDoc description
public class ClusterCleanerUtilities {

    private List<ArrayList<Hit>> sortedHits;

    public ClusterCleanerUtilities() {
        List<ArrayList<Hit>> sortdHits = new ArrayList<ArrayList<Hit>>();
        for(int l = 0; l < 6; l++) sortdHits.add(new ArrayList<Hit>());
        sortedHits = sortdHits;
    }

    /**
     * Pattern Recognition step for identifying clusters in a clump:
     * Find the points that are  consistent with belonging to the same cluster. This step precedes
     * the initial estimates of the track segments which require further refining. The method
     * employed is that of Hough Transforms. Define the dimension of the r-theta accumulator array
     * used for pattern recognition of (rho, phi) points.
     *
     * @param clus              the fitted cluster. This cluster is examined for overlaps and tracks
     * @param nextClsStartIndex the index of the next cluster in the splitted cluster
     * @param cf                instance of the cluster fitter class
     * @return                  a list of fitted clusters
     */
    public List<FittedCluster> clusterSplitter(FittedCluster clus,
                                               int nextClsStartIndex,
                                               ClusterFitter cf) {
        /*
        The principle of Hough Transform in pattern recognition is as follows:

        For every point (rho, phi) on a line there exists an infinite number of lines that go
        through it.  Each such line can be parametrized with parameters (r, theta) such that
        r = rho * cos(theta) + phi * sin(theta), where theta is the polar angle of the line which is
        perpendicular to it and intersects the origin, and r is the distance between that line and
        the origin, for a given theta. Hence a point in (rho, phi) parameter space corresponds to a
        sinusoidal curve in (r, theta) parameter space, which is the so-called Hough-transform
        space. Points that fall on a line in (rho, phi) space correspond to sinusoidal curves that
        intersect at a common point in Hough-transform space. Remapping this point in (rho, phi)
        space yields the line that contains the points in the original line.

        This method is a pattern recognition tool used to select groups of points belonging to the
        same shape, such as a line or a circle. To find the ensemble of points belonging to a line
        in the original space the Hough transform algorithm makes use of an array, called an
        accumulator.

        The principle of the accumulator is a counting method. The dimensions of the array is equal
        to the number of parameters in Hough transform space, which is 2 corresponding to the
        (r, theta) pair in our particular case.

        The bin size in the array are finite intervals in r and theta, which are called accumulator
        cells. The bin content of cells along the curve of discretized (r, theta) values get
        incremented. The cell with the highest count corresponds to the intersection of the curves.

        This is a numerical method to find the intersection of any number of curves. Once the
        accumulator array has been filled with all (r, theta) points peaks and their associated
        (rho, phi) points are determined. From these, sets of points belonging to common lines can be
        determined. This is a preliminary pattern recognition method used to identify reconstructed
        hits belonging to the same track-segment.
        */
        int N_t = 180;

        // From this calculate the bin size in the theta accumulator array
        double ThetaMin = 0.;
        double ThetaMax = 2. * Math.PI;
        double SizeThetaBin = (ThetaMax - ThetaMin) / ((double) N_t);

        // Define the dimension of the r accumulator array
        int N_r = 130;

        // From this calculate the bin size in the theta accumulator array
        double RMin = -130;
        double RMax = 130;

        int[][] R_Phi_Accumul;
        R_Phi_Accumul = new int[N_r][N_t];

        // Cache the cos and sin theta values [for performance improvement]
        double[] cosTheta_RPhi_array;
        double[] sinTheta_RPhi_array;

        // The values corresponding to the peaks in the array
        double[] binrMaxR_Phi;
        double[] bintMaxR_Phi;
        binrMaxR_Phi = new double[N_r * N_t];
        bintMaxR_Phi = new double[N_r * N_t];

        cosTheta_RPhi_array = new double[N_t];
        sinTheta_RPhi_array = new double[N_t];

        for (int j_t = 0; j_t < N_t; j_t++) {
            // theta_j in the middle of the bin :
            double theta_j = ThetaMin + (0.5 + j_t) * SizeThetaBin;
            cosTheta_RPhi_array[j_t] = Math.cos(theta_j);
            sinTheta_RPhi_array[j_t] = Math.sin(theta_j);
        }

        // Loop over points to fill the accumulator arrays
        for (int i = 0; i < clus.size(); i++) {

            double rho = clus.get(i).get_lX();
            double phi = clus.get(i).get_lY();

            // Fill the accumulator arrays
            for (int j_t = 0; j_t < N_t; j_t++) {
                // r_j corresponding to theta_j
                double r_j = rho * cosTheta_RPhi_array[j_t] + phi * sinTheta_RPhi_array[j_t];

                // This value of r_j falls into the following bin in the r array:
                int j_r = (int) Math.floor(N_r * (r_j - RMin) / (float) (RMax - RMin));

                // Increase this accumulator cell:
                R_Phi_Accumul[j_r][j_t]++;
            }
        }

        // Loop over accumulator array to find peaks (allows for more than one peak for multiple
        //     tracks).
        // The accumulator cell count must be at least half the total number of hits
        // Make binrMax, bintMax arrays to allow for more than one peak
        int threshold = Constants.DC_MIN_NLAYERS;
        int nbPeaksR_Phi = 0;

        // 1st find the peaks in the R_Phi accumulator array
        for (int ibinr1 = 0; ibinr1 < N_r; ibinr1++) {
            for (int ibint1 = 0; ibint1 < N_t; ibint1++) {

                // Find the peak
                if (R_Phi_Accumul[ibinr1][ibint1] >= Constants.DC_MIN_NLAYERS) {
                    if (R_Phi_Accumul[ibinr1][ibint1] > threshold) {
                        threshold = R_Phi_Accumul[ibinr1][ibint1];
                    }

                    binrMaxR_Phi[nbPeaksR_Phi] = ibinr1;
                    bintMaxR_Phi[nbPeaksR_Phi] = ibint1;
                    nbPeaksR_Phi++;
                }
            }
        }

        // For a given Maximum value of the accumulator, find the set of points associated with it.
        // For this, loop again over all the points.
        List<FittedCluster> splitclusters = new ArrayList<FittedCluster>();

        for (int p = nbPeaksR_Phi - 1; p > -1; p--) {

            // Make a new cluster
            FittedCluster newClus = new FittedCluster(clus.getBaseCluster());

            for (int i = 0; i < clus.size(); i++) {
                double rho = clus.get(i).get_X();
                double phi = clus.get(i).get_lY();

                for (int j_t = 0; j_t < N_t; j_t++) {
                    // r_j corresponding to theta_j:
                    double r_j = rho * cosTheta_RPhi_array[j_t] + phi * sinTheta_RPhi_array[j_t];
                    // This value of r_j falls into the following bin in the r array:
                    int j_r = (int) Math.floor(N_r * (r_j - RMin) / (float) (RMax - RMin));

                    // Match bins:
                    if (j_r == binrMaxR_Phi[p] && j_t == bintMaxR_Phi[p]) newClus.add(clus.get(i));
                }
            }

            List<Hit> contigArrayOfHits = new ArrayList<Hit>(); // Contiguous cluster

            boolean passCluster = true;
            for (int l = 1; l <= Constants.NLAYR; l++) {
                for (int i = 0; i < newClus.size(); i++) {
                    if (newClus.get(i).get_Layer() == l) {
                        contigArrayOfHits.add(newClus.get(i));
                    }
                }
            }
            for (int i = 0; i < contigArrayOfHits.size() - 1; i++) {
                // If there is a gap, do not include in list
                if (contigArrayOfHits.get(i + 1).get_Layer()
                        - contigArrayOfHits.get(i).get_Layer() > 1) {
                    passCluster = false;
                }
            }
            // Require 4 layers to make a cluster
            if (countLayersInCluster(contigArrayOfHits) < Constants.DC_MIN_NLAYERS) {
                passCluster = false;
            }

            // Require consistency with line
            cf.SetFitArray(newClus, "LC");
            cf.Fit(newClus, true);

            if (newClus.get_fitProb() < 0.9) passCluster = false;
            if (!(splitclusters.contains(newClus)) && passCluster) splitclusters.add(newClus);
        }

        // Make new clusters
        List<FittedCluster> selectedClusList = new ArrayList<FittedCluster>();

        int newcid = nextClsStartIndex;
        for (FittedCluster cluster : splitclusters) {
            cluster.set_Id(newcid++);
            cf.SetFitArray(cluster, "LC");
            cf.Fit(cluster, true);

            FittedCluster bestCls = overlappingClusterResolver(cluster, splitclusters);

            if (bestCls != null) {
                if (!(selectedClusList.contains(bestCls))) selectedClusList.add(bestCls);
            }
        }

        int splitclusId = 1;
        if (selectedClusList.size() != 0) {
            for (FittedCluster cl : selectedClusList) {
                cl.set_Id(clus.get_Id() * 1000 + splitclusId);
                splitclusId++;
            }
        }

        if (selectedClusList.size() == 0) {
            // If the splitting fails, then return the original cluster
            selectedClusList.add(clus);
        }
        return selectedClusList;
    }

    /**
     * Sorts a list of hits by layer.
     * @param DCHits   list of hits to be sorted
     * @param sector   NOTE: lacks description
     * @param superlyr NOTE: lacks description
     * @return         sorted list of hits
     */
    public List<List<Hit>> sortByLayerList(List<Hit> DCHits, int sector, int superlyr) {

        List<List<Hit>> hitsinlayr_array = new ArrayList<List<Hit>>();
        int nlayr = 6;
        for (int l = 0; l < nlayr; l++) {
            List<Hit> hitsinlayr = new ArrayList<Hit>();
            for (Hit hitInList : DCHits) {
                if (hitInList != null) {
                    if (hitInList.get_Layer() == l + 1
                            && hitInList.get_Sector() == sector
                            && hitInList.get_Superlayer() == superlyr) {
                        hitsinlayr.add(hitInList);
                    }
                }
            }
            hitsinlayr_array.add(l, hitsinlayr);
        }
        return hitsinlayr_array;
    }

    /**
     * Finds and returns the number of layers hit at a certain wire coordinate.
     * @param hits_inlayer the hits in a given layer
     * @return             number of layers hit
     */
    public int countLayersHit(Hit[] hits_inlayer) {
        int nlayr = 6;
        Hit[] allhits_inlayer = new Hit[nlayr];
        allhits_inlayer = hits_inlayer;

        int nlayers_hit = 0;
        for (int la = 0; la < nlayr; la++) {
            if (allhits_inlayer[la] != null) nlayers_hit++;
        }
        return nlayers_hit;
    }

    /**
     * Counts the layers in a cluster.
     * @param hitsInClus the hits in a cluster
     * @return           the number of layers in a cluster
     */
    public int countLayersInCluster(List<Hit> hitsInClus) {
        // Count hits in each layer
        int nlayr = 6;
        int[] nlayers = new int[nlayr];
        for (int l = 0; l < nlayr; l++) {
            nlayers[l] = 0;
            for (int h = 0; h < hitsInClus.size(); h++) {
                if (hitsInClus.get(h).get_Layer() == l + 1) nlayers[l]++;
            }
        }

        // Count number of layers hit
        int nlayers_hit = 0;
        for (int l = 0; l < nlayr; l++) {
            if (nlayers[l] > 0) nlayers_hit++;
        }

        return nlayers_hit;
    }

    /**
     * Resolves the Left-Right ambiguity in a cluster.
     * @param fClus      the cluster
     * @param cf         a cluster fitter instance
     * @param tab        NOTE: parameter lacks description
     * @param DcDetector DC Detector Geometry
     * @param tde        NOTE: parameter lacks description
     * @return           NOTE: return data lacks description
     */
    public FittedCluster LRAmbiguityResolver(FittedCluster fClus,
                                             ClusterFitter cf,
                                             IndexedTable tab,
                                             DCGeant4Factory DcDetector,
                                             TimeToDistanceEstimator tde) {
        int index = 0;
        for(FittedHit hit: fClus) {
            if (hit.get_Doca() < 0.4 * hit.get_CellSize()) hit.set_LeftRightAmb(0);
            if (hit.get_LeftRightAmb() == 0) index++;
        }
        if (index == 0) return fClus; // Cluster OK
        if (index > 6)  return null;  // Unresolveable cluster...
        int arraySize = (int) Math.pow(2, (double) index);
        ArrayList<FittedCluster> arrayOfClus = new ArrayList<FittedCluster>(arraySize);

        // Pass all acceptable clusters
        FittedCluster okClus = new FittedCluster(fClus.getBaseCluster());
        for (FittedHit hit : fClus) {
            if (hit.get_LeftRightAmb() != 0) okClus.add(hit);
        }

        // Filter all other clusters
        FittedCluster notLRClus = new FittedCluster(fClus.getBaseCluster());
        for (FittedHit hit : fClus) {
            if (hit.get_LeftRightAmb() == 0) notLRClus.add(hit);
        }

        // Make combinatorials
        FittedCluster tnLRc = new FittedCluster(fClus.getBaseCluster());        // totNotLRClus
        FittedCluster pnLRc = new FittedCluster(fClus.getBaseCluster());        // posNotLRClus
        FittedCluster nnLRc = new FittedCluster(fClus.getBaseCluster());        // negNotLRClus

        for (FittedHit hit : notLRClus) {

            FittedHit newhitPos = new FittedHit(hit.get_Sector(), hit.get_Superlayer(),
                                                hit.get_Layer(),  hit.get_Wire(),
                                                hit.get_TDC(),    hit.get_Id());
            newhitPos.set_Doca(hit.get_Doca());
            newhitPos.set_DocaErr(hit.get_DocaErr());
            newhitPos.setT0(hit.getT0());
            newhitPos.set_Beta(hit.get_Beta());
            newhitPos.setB(hit.getB());
            newhitPos.set_DeltaTimeBeta(hit.get_DeltaTimeBeta());
            newhitPos.setTStart(hit.getTStart());
            newhitPos.setTProp(hit.getTProp());
            newhitPos.setTFlight(hit.getTFlight());
            newhitPos.set_Time(hit.get_Time());
            newhitPos.set_Id(hit.get_Id());
            newhitPos.set_TrkgStatus(0);
            newhitPos.calc_CellSize(DcDetector);
            newhitPos.set_LeftRightAmb(1);
            // Assume the track angle is parallel to the layer, so that cosTrkAng = 1
            newhitPos.updateHitPositionWithTime(1, hit.getB(), tab, DcDetector, tde);

            newhitPos.set_AssociatedClusterID(hit.get_AssociatedClusterID());
            newhitPos.set_AssociatedHBTrackID(hit.get_AssociatedHBTrackID());

            FittedHit newhitNeg = new FittedHit(hit.get_Sector(), hit.get_Superlayer(),
                                                hit.get_Layer(),  hit.get_Wire(),
                                                hit.get_TDC(),    hit.get_Id());
            newhitNeg.set_Doca(hit.get_Doca());
            newhitNeg.set_DocaErr(hit.get_DocaErr());
            newhitNeg.setT0(hit.getT0());
            newhitNeg.set_Beta(hit.get_Beta());
            newhitNeg.setB(hit.getB());
            newhitNeg.set_DeltaTimeBeta(hit.get_DeltaTimeBeta());
            newhitNeg.setTStart(hit.getTStart());
            newhitNeg.setTProp(hit.getTProp());
            newhitNeg.setTFlight(hit.getTFlight());
            newhitNeg.set_Time(hit.get_Time());
            newhitNeg.set_Id(hit.get_Id());
            newhitNeg.set_TrkgStatus(0);
            newhitNeg.calc_CellSize(DcDetector);
            newhitNeg.set_LeftRightAmb(-1);
            // Assume the track angle is parallel to the layer
            newhitNeg.updateHitPositionWithTime(1, hit.getB(), tab, DcDetector, tde);

            newhitNeg.set_AssociatedClusterID(hit.get_AssociatedClusterID());
            newhitNeg.set_AssociatedHBTrackID(hit.get_AssociatedHBTrackID());

            tnLRc.add(newhitNeg);
            tnLRc.add(newhitPos);

            pnLRc.add(newhitPos);
            nnLRc.add(newhitNeg);
        }

        Collections.sort(tnLRc);

        if (index == 1) {
            arrayOfClus.add(pnLRc);
            arrayOfClus.add(nnLRc);
            arrayOfClus.get(0).addAll(okClus);
            arrayOfClus.get(1).addAll(okClus);
        }
        else if (index == 2) {
            for (int i1 = 0; i1 < tnLRc.size(); i1++) {
                for (int i2 = 2; i2 < tnLRc.size(); i2++) {
                    if (tnLRc.get(i1).get_Id() == tnLRc.get(i2).get_Id()) {
                        continue;
                    }
                    FittedCluster newClus = new FittedCluster(fClus.getBaseCluster());
                    newClus.addAll(okClus);
                    newClus.add(tnLRc.get(i1));
                    newClus.add(tnLRc.get(i2));
                    arrayOfClus.add(newClus);
                }
            }
        }

        else if (index == 3) {
            for (int i1 = 0; i1 < tnLRc.size(); i1++) {
                for (int i2 = 2; i2 < tnLRc.size(); i2++) {
                    for (int i3 = 4; i3 < tnLRc.size(); i3++) {
                        if ((tnLRc.get(i1).get_Id() == tnLRc.get(i2).get_Id())
                                || (tnLRc.get(i1).get_Id() == tnLRc.get(i3).get_Id())
                                || (tnLRc.get(i2).get_Id() == tnLRc.get(i3).get_Id())) {
                            continue;
                        }
                        FittedCluster newClus = new FittedCluster(fClus.getBaseCluster());
                        newClus.addAll(okClus);
                        newClus.add(tnLRc.get(i1));
                        newClus.add(tnLRc.get(i2));
                        newClus.add(tnLRc.get(i3));
                        arrayOfClus.add(newClus);
                    }
                }
            }
        }

        else if (index == 4) {
            for (int i1 = 0; i1 < tnLRc.size(); i1++) {
                for (int i2 = 2; i2 < tnLRc.size(); i2++) {
                    for (int i3 = 4; i3 < tnLRc.size(); i3++) {
                        for (int i4 = 6; i4 < tnLRc.size(); i4++) {
                            if ((tnLRc.get(i1).get_Id() == tnLRc.get(i2).get_Id())
                                    || (tnLRc.get(i1).get_Id() == tnLRc.get(i3).get_Id())
                                    || (tnLRc.get(i1).get_Id() == tnLRc.get(i4).get_Id())
                                    || (tnLRc.get(i2).get_Id() == tnLRc.get(i3).get_Id())
                                    || (tnLRc.get(i2).get_Id() == tnLRc.get(i4).get_Id())
                                    || (tnLRc.get(i3).get_Id() == tnLRc.get(i4).get_Id())) {
                                continue;
                            }
                            FittedCluster newClus = new FittedCluster(fClus.getBaseCluster());
                            newClus.addAll(okClus);
                            newClus.add(tnLRc.get(i1));
                            newClus.add(tnLRc.get(i2));
                            newClus.add(tnLRc.get(i3));
                            newClus.add(tnLRc.get(i4));
                            arrayOfClus.add(newClus);
                        }
                    }
                }
            }
        }

        else if (index == 5) {
            for (int i1 = 0; i1 < tnLRc.size(); i1++) {
                for (int i2 = 2; i2 < tnLRc.size(); i2++) {
                    for (int i3 = 4; i3 < tnLRc.size(); i3++) {
                        for (int i4 = 6; i4 < tnLRc.size(); i4++) {
                            for (int i5 = 8; i5 < tnLRc.size(); i5++) {
                                if ((tnLRc.get(i1).get_Id() == tnLRc.get(i2).get_Id())
                                        || (tnLRc.get(i1).get_Id() == tnLRc.get(i3).get_Id())
                                        || (tnLRc.get(i1).get_Id() == tnLRc.get(i4).get_Id())
                                        || (tnLRc.get(i1).get_Id() == tnLRc.get(i5).get_Id())
                                        || (tnLRc.get(i2).get_Id() == tnLRc.get(i3).get_Id())
                                        || (tnLRc.get(i2).get_Id() == tnLRc.get(i4).get_Id())
                                        || (tnLRc.get(i2).get_Id() == tnLRc.get(i5).get_Id())
                                        || (tnLRc.get(i3).get_Id() == tnLRc.get(i4).get_Id())
                                        || (tnLRc.get(i3).get_Id() == tnLRc.get(i5).get_Id())
                                        || (tnLRc.get(i4).get_Id() == tnLRc.get(i5).get_Id())) {
                                    continue;
                                }
                                FittedCluster newClus = new FittedCluster(fClus.getBaseCluster());
                                newClus.addAll(okClus);
                                newClus.add(tnLRc.get(i1));
                                newClus.add(tnLRc.get(i2));
                                newClus.add(tnLRc.get(i3));
                                newClus.add(tnLRc.get(i4));
                                newClus.add(tnLRc.get(i5));
                                arrayOfClus.add(newClus);
                            }
                        }
                    }
                }
            }
        }

        else if (index == 6) {
            for (int i1 = 0; i1 < tnLRc.size(); i1++) {
                for (int i2 = 2; i2 < tnLRc.size(); i2++) {
                    for (int i3 = 4; i3 < tnLRc.size(); i3++) {
                        for (int i4 = 6; i4 < tnLRc.size(); i4++) {
                            for (int i5 = 8; i5 < tnLRc.size(); i5++) {
                                for (int i6 = 10; i6 < tnLRc.size(); i6++) {
                                    if ((tnLRc.get(i1).get_Id() == tnLRc.get(i2).get_Id())
                                            || (tnLRc.get(i1).get_Id() == tnLRc.get(i3).get_Id())
                                            || (tnLRc.get(i1).get_Id() == tnLRc.get(i4).get_Id())
                                            || (tnLRc.get(i1).get_Id() == tnLRc.get(i5).get_Id())
                                            || (tnLRc.get(i1).get_Id() == tnLRc.get(i6).get_Id())
                                            || (tnLRc.get(i2).get_Id() == tnLRc.get(i3).get_Id())
                                            || (tnLRc.get(i2).get_Id() == tnLRc.get(i4).get_Id())
                                            || (tnLRc.get(i2).get_Id() == tnLRc.get(i5).get_Id())
                                            || (tnLRc.get(i2).get_Id() == tnLRc.get(i6).get_Id())
                                            || (tnLRc.get(i3).get_Id() == tnLRc.get(i4).get_Id())
                                            || (tnLRc.get(i3).get_Id() == tnLRc.get(i5).get_Id())
                                            || (tnLRc.get(i3).get_Id() == tnLRc.get(i6).get_Id())
                                            || (tnLRc.get(i4).get_Id() == tnLRc.get(i5).get_Id())
                                            || (tnLRc.get(i4).get_Id() == tnLRc.get(i6).get_Id())
                                            || (tnLRc.get(i5).get_Id() == tnLRc.get(i6).get_Id())) {
                                        continue;
                                    }
                                    FittedCluster newClus = new FittedCluster(fClus.getBaseCluster());
                                    newClus.addAll(okClus);
                                    newClus.add(tnLRc.get(i1));
                                    newClus.add(tnLRc.get(i2));
                                    newClus.add(tnLRc.get(i3));
                                    newClus.add(tnLRc.get(i4));
                                    newClus.add(tnLRc.get(i5));
                                    newClus.add(tnLRc.get(i6));
                                    arrayOfClus.add(newClus);
                                }
                            }
                        }
                    }
                }
            }
        }

        return cf.BestClusterSelector(arrayOfClus, "TSC");
    }

    /**
     * NOTE: method lacks description
     * @param Clus       the cluster
     * @param cf         a cluster fitter instance
     * @param tab        NOTE: parameter lacks description
     * @param DcDetector DC Detector Geometry
     * @param tde        NOTE: parameter lacks description
     * @return           NOTE: return data lacks description
     */
    public FittedCluster secondariesRemover(FittedCluster clus,
                                            ClusterFitter cf,
                                            IndexedTable tab,
                                            DCGeant4Factory DcDetector,
                                            TimeToDistanceEstimator tde) {
        Collections.sort(clus);

        ArrayList<ArrayList<FittedHit>> sortedHits = new ArrayList<ArrayList<FittedHit>>(6);
        // Initialize
        for (int i = 0; i < 6; i++) sortedHits.add(new ArrayList<FittedHit>());

        for (int i = 0; i < clus.size(); i++) {
            FittedHit hit = clus.get(i);
            sortedHits.get(hit.get_Layer() - 1).add(hit);
        }

        ArrayList<FittedCluster> clusters = new ArrayList<FittedCluster>();

        ArrayList<ArrayList<FittedHit>> hitsInSameLayerLists = new ArrayList<ArrayList<FittedHit>>();
        ArrayList<ArrayList<FittedHit>> hitsInClusCandLists = new ArrayList<ArrayList<FittedHit>>();
        ArrayList<FittedHit> baseClusterHits = new ArrayList<FittedHit>();

        for (int i = 0; i < 6; i++) {
            ArrayList<FittedHit> hitsInLayer = sortedHits.get(i);
            if (hitsInLayer.size() == 0) continue;
            if (hitsInLayer.size() == 1) {
                // Save all good hits to base cluster
                baseClusterHits.addAll(hitsInLayer);
                for (int j = 0; j < hitsInLayer.size(); j++) {
                    hitsInLayer.get(j).set_LeftRightAmb(0);
                    hitsInLayer.get(j).updateHitPositionWithTime(1, hitsInLayer.get(j).getB(),
                                                                 tab, DcDetector, tde);
                }
            }
            if (hitsInLayer.size() == 2) {
                double docaSum = 0;
                for (int j = 0; j < hitsInLayer.size(); j++) {
                    docaSum += hitsInLayer.get(j).get_Doca();
                }

                double passingCut = 1.5 * hitsInLayer.get(0).get_CellSize() ;

                double hit1doca = 0;
                double hit2doca = 0;
                if (hitsInLayer.get(0).get_Doca() > hitsInLayer.get(1).get_Doca()) {
                    hit1doca = hitsInLayer.get(0).get_Doca();
                    hit2doca = hitsInLayer.get(1).get_Doca();
                } else {
                    hit1doca = hitsInLayer.get(1).get_Doca();
                    hit2doca = hitsInLayer.get(0).get_Doca();
                }
                double passingCut2 = 0.75;
                if (hit2doca/hit1doca < passingCut2
                        || (hit2doca/hit1doca > passingCut2 && docaSum < passingCut)) {
                    // reset LR to 0
                    for (int j = 0; j < hitsInLayer.size(); j++) {
                        hitsInLayer.get(j).set_LeftRightAmb(0);
                        hitsInLayer.get(j).updateHitPositionWithTime(1, hitsInLayer.get(j).getB(),
                                                                     tab, DcDetector, tde);
                    }
                    hitsInSameLayerLists.add(hitsInLayer);
                } else {
                    // Save all good hits to base cluster
                    baseClusterHits.addAll(hitsInLayer);
                }
            }
        }

        int nbLyr = hitsInSameLayerLists.size();
        if      (nbLyr == 0) return clus;
        else if (nbLyr > 0) {
            for (int clusIdx = 0; clusIdx < Constants.CombArray.get(nbLyr - 1).length; clusIdx++) {
                ArrayList<FittedHit> hitsInClusterCand = new ArrayList<FittedHit>();
                hitsInClusterCand.addAll(baseClusterHits);

                for (int k = 0; k < Constants.CombArray.get(nbLyr - 1)[clusIdx].length; k++) {
                    hitsInClusterCand.add(hitsInSameLayerLists.get(k).get(
                            Constants.CombArray.get(nbLyr - 1)[clusIdx][k]));
                }
                hitsInClusCandLists.add(hitsInClusterCand);
            }
        }
        for (int i = 0; i < hitsInClusCandLists.size(); i++) {
            FittedCluster newClus = new FittedCluster(clus.getBaseCluster());
            for (int i1 = 0; i1 < newClus.size(); i1++) newClus.remove(i1);

            newClus.addAll(hitsInClusCandLists.get(i));
            clusters.add(newClus);
        }

        FittedCluster BestCluster = cf.BestClusterSelector(clusters, "LC");

        return BestCluster;
    }

    /**
     * A method to select the largest cluster among a set of clusters with 4 or more of overlaping
     * hits.
     * @param thisclus the cluster to be compared to a list of other clusters
     * @param clusters the list of clusters
     * @return         the selected cluster
     */
    public FittedCluster overlappingClusterResolver(FittedCluster thisclus,
                                                    List<FittedCluster> clusters) {

        List<FittedCluster> overlapingClusters = new ArrayList<FittedCluster>();

        for (FittedCluster cls : clusters) {
            List<FittedHit> hitOvrl = new ArrayList<FittedHit>();
            for (FittedHit hit : thisclus) {
                if (cls.contains(hit)) {
                    if (!(hitOvrl.contains(hit))) hitOvrl.add(hit);
                }
            }

            // Test
            boolean passCls = true;
            for (FittedCluster ovr : overlapingClusters) {
                if (ovr.get_Id() == cls.get_Id()) {
                    passCls = false;
                    continue;
                }

                // Ensure that the lines are consistent
                if (Math.abs(ovr.get_clusterLineFitSlope() - cls.get_clusterLineFitSlope()) > 0.2) {
                    passCls = false;
                }
            }
            if (hitOvrl.size() < 3) passCls = false;
            if (passCls) overlapingClusters.add(cls);

        }
        Collections.sort(overlapingClusters);

        // Return the largest cluster.
        return overlapingClusters.get(0);
    }

    /**
     * Prunes the input hit list to remove noise candidates; the algorithm finds contiguous hits in
     * a layer (column) and removes hits according to the number (Nc) of such contiguous hits in a
     * given layer.
     * If Nc = 3, keep only the middle hit.
     * If Nc = 4, keep only the first and last hit in that column.
     * if Nc > 4, keep the first 2 and last 2 hits in that column.
     * if Nc > 10 remove all hits in that column.
     * @param hits the unfitted hits
     * @return     the pruned list of hits
     */
    public List<Hit> hitListPruner(List<Hit> hits) {
        for (int l = 0; l < 6; l++) sortedHits.get(l).clear();
        for (int i = 0; i < hits.size() ; i++) {
            sortedHits.get(hits.get(i).get_Layer()-1).add(hits.get(i));
        }

        for (int l = 0; l < 6; l++) {
            if (sortedHits.get(l).size() == 10) {
                sortedHits.get(l).removeAll(sortedHits.get(l));
            }
            else if (sortedHits.get(l).size() > 4) {
                ArrayList<Hit> rmHits = (ArrayList<Hit>) sortedHits.get(l).clone();
                ArrayList<Hit> kHits = new ArrayList<Hit>();
                kHits.add(sortedHits.get(l).get(0));
                kHits.add(sortedHits.get(l).get(1));
                kHits.add(sortedHits.get(l).get(sortedHits.get(l).size()-1));
                kHits.add(sortedHits.get(l).get(sortedHits.get(l).size()-2));
                rmHits.removeAll(kHits);
                sortedHits.get(l).removeAll(rmHits);
            }
            else if (sortedHits.get(l).size() > 2) {
                ArrayList<Hit> rmHits = (ArrayList<Hit>) sortedHits.get(l).clone();
                ArrayList<Hit> kHits  = new ArrayList<Hit>();
                kHits.add(sortedHits.get(l).get(0));
                kHits.add(sortedHits.get(l).get(sortedHits.get(l).size()-1));
                rmHits.removeAll(kHits);
                sortedHits.get(l).removeAll(rmHits);
            }
        }

        hits.clear();
        for(int l = 0; l < 6; l++) {
            if (sortedHits.get(l).size() > 0) hits.addAll(sortedHits.get(l));
        }
        return hits;
    }

    /**
     * Prunes the isolated hits in a cluster.
     * @param clus un-prunned cluster
     * @return     new contiguous cluster
     */
    public FittedCluster isolatedHitsPruner(FittedCluster clus) {

        int min = 1000;
        int max = -1000;
        for (int i = 0; i < clus.size(); i++) {
            if (clus.get(i).get_Wire() <= min) min = clus.get(i).get_Wire();
            if (clus.get(i).get_Wire() >= max) max = clus.get(i).get_Wire();
        }
        min -= 1;
        max += 1;
        int wireRange = max - min + 1;
        FittedHit[][] HitArray = new FittedHit[6][wireRange];

        for (int i = 0; i < clus.size(); i++) {
            int wi = clus.get(i).get_Wire() - min;
            if (wi < -1) continue;

            int la = clus.get(i).get_Layer() - 1;
            if (wi >= 0 && wi < wireRange) HitArray[la][wi] = clus.get(i);
        }

        Cluster newCluster = new Cluster(clus.get_Sector(), clus.get_Superlayer(), clus.get_Id());
        FittedCluster fcluster = new FittedCluster(newCluster);

        for (int i = 0; i < clus.size(); i++) {
            int wire = clus.get(i).get_Wire() - min + 1;
            int layer = clus.get(i).get_Layer();
            if (layer == 1) {
                // Look for neighbor in next layer
                if (HitArray[layer][wire - 1]            != null
                        || HitArray[layer][wire - 2]     != null
                        || HitArray[layer][wire]         != null
                        || HitArray[layer - 1][wire - 2] != null
                        || HitArray[layer - 1][wire]     != null) {
                    fcluster.add(clus.get(i));
                }
            }
            else if (layer == 6) {
            // Look for neighbor in previous layer
                if (HitArray[layer - 2][wire - 1]        != null
                        || HitArray[layer - 2][wire - 2] != null
                        || HitArray[layer - 2][wire]     != null
                        || HitArray[layer - 1][wire - 2] != null
                        || HitArray[layer - 1][wire]     != null) {
                    fcluster.add(clus.get(i));
                }
            }
            else {
                // Look for neighbor in next and previous layers
                if (HitArray[layer][wire - 1]            != null
                        || HitArray[layer][wire - 2]     != null
                        || HitArray[layer][wire]         != null
                        || HitArray[layer - 2][wire - 1] != null
                        || HitArray[layer - 2][wire - 2] != null
                        || HitArray[layer - 2][wire]     != null
                        || HitArray[layer - 1][wire - 2] != null
                        || HitArray[layer - 1][wire]     != null) {
                    fcluster.add(clus.get(i));
                }
            }
        }

        return fcluster;
    }

    /**
     * Removes out of time hits.
     * @param fclus     fitted cluster to be filtered
     * @param removeHit flag telling if hits should be removed
     */
    public void outOfTimersRemover(FittedCluster fClus, boolean removeHit) {
        for (int i = 0; i < fClus.size(); i++) {
            if (fClus.get(i).get_OutOfTimeFlag() == true) {
                if (removeHit == true) fClus.remove(i);
                else fClus.get(i).set_Doca(fClus.get(i).get_CellSize());
            }
        }
    }
}
