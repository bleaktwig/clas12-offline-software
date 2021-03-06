package org.jlab.rec.dc.cluster;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.CopyOnWriteArrayList;

import org.jlab.detector.geant4.v2.DCGeant4Factory;
import org.jlab.io.evio.EvioDataBank;
import org.jlab.io.evio.EvioDataEvent;
import org.jlab.utils.groups.IndexedTable;
import org.jlab.rec.dc.Constants;
import org.jlab.rec.dc.hit.Hit;
import org.jlab.rec.dc.hit.FittedHit;
import org.jlab.rec.dc.timetodistance.TimeToDistanceEstimator;

/**
 * A hit pruning algorithm to reject noise that gives a pattern of hits that are continguous in the
 * same layer. The algorithm first puts the hits in arrays according to their layer and wire number.
 * Each such array contains all the hits in the same layer. The algorithm then collects groups of
 * contiguous hits into a list of hits. The n-first and n-last hits in the list are kept, with all
 * other hits inbetween pruned. The value of n depends on the size of the list. A loose clustering
 * algorithm loops over all superlayers, in a sector and finds groups of hits with contiguous wire
 * index numbers. These clusters (called clumps of hits) are delimited by layers with no hits at a
 * particular wire coordinate, and are then refined using fits to their respective wire indexes as a
 * function of layer number to identify parallel tracks or overlapping track candidates.
 *
 * @author ziegler
 */
public class ClusterFinder {

    public ClusterFinder() {}

    // the loop is done over sector and superlayers
    // idx        = superlayer*sector + superlayer
    // sector     = idx/nsect + 1 (starting at 1)
    // superlayer = idx%nsect + 1 (starting at 1)
    int nsect = Constants.NSECT;
    int nslay = Constants.NSLAY;
    int nlayr = Constants.NLAYR;
    int nwire = Constants.NWIRE;

    // 3-dimentional array of hits as an array:
    // [total_nb_sectors*total_nb_superlayers][total_nb_wires][total_nb_layers]
    private Hit[][][] HitArray = new Hit[nsect * nslay][nwire][nlayr];

    public Hit[][][] getHitArray() {return HitArray;}
    public void setHitArray(Hit[][][] hitArray) {HitArray = hitArray;}

    /**
     * Fills 3-dimentional array of hits from input hits
     * @param hits              the unfitted hit
     * @param rejectLayer
     */
    private void fillHitArray(List<Hit> hits, int rejectLayer) {

        Hit[][][] hitArray = new Hit[nsect * nslay][nwire][nlayr];

        // initializing non-zero Hit Array entries with valid hits
        for (Hit hit : hits) {
            if (passHitSelection(hit) && hit.get_Layer() != rejectLayer) {
                int ssl = (hit.get_Sector() - 1) * nsect + (hit.get_Superlayer() - 1);
                int wi = hit.get_Wire() - 1;
                int la = hit.get_Layer() - 1;

                if (wi >= 0 && wi < nwire) hitArray[ssl][wi][la] = hit;
            }
        }
        this.setHitArray(hitArray);
    }

    /**
     * Finds clumps. A clump is a cluster that is not filtered for noise.
     * @param allhits the list of unfitted hits
     * @param ct
     * @return List of clusters
     */
    public List<Cluster> findClumps(List<Hit> allhits, ClusterCleanerUtilities ct) {
        Collections.sort(allhits);

        List<Cluster> clumps = new ArrayList<Cluster>();

        // Looping over each superlayer in each sector
        // Each superlayer is treated independently
        int cid = 1; // cluster id, will increment with each new good cluster

        for (int ssl = 0; ssl < nsect * nslay; ssl++) {
            // For each ssl, a loop over the wires is done to define clusters.
            // Clusters are delimited by layers with no hits at a particular wire coordinate.

            int wi = 0; // wire number in the loop
            // Looping over all physical wires
            while (wi < nwire) {
                // if there's a hit in at least one layer, it's a cluster candidate
                if (ct.countLayersHit(HitArray[ssl][wi]) != 0) {
                    List<Hit> hits = new ArrayList<Hit>();

                    // Adding all hits in this and all the subsequent
                    // Wires until there's a wire with no layers hit
                    while (ct.countLayersHit(HitArray[ssl][wi]) > 0 && wi < nwire) {
                        // Looping over all physical wires
                        for (int la = 0; la < nlayr; la++) {
                            if (HitArray[ssl][wi][la] != null) hits.add(HitArray[ssl][wi][la]);
                        }
                        wi++;
                    }

                    // Need at least MIN_NLAYERS
                    if (ct.countLayersInCluster(hits) >= Constants.DC_MIN_NLAYERS) {
                        // cluster constructor DCCluster(hit.sector,hit.superlayer, cid)
                        Cluster this_cluster = new Cluster((int) (ssl / nsect) + 1,
                                                           (int) (ssl % nsect) + 1,
                                                           cid++);
                        this_cluster.addAll(hits);
                        clumps.add(this_cluster);
                    }
                }
                wi++; // if no hits are added, check for the next wire coordinate
            }
        }
        return clumps;
    }

    /**
     * Hit-based tracking linear fits to the wires is done to determine the cluster, resulting in a
     * fitted cluster.
     * @param allhits    the list of unfitted hits
     * @param ct         ClusterCleanerUtilities instance
     * @param cf         ClusterFitter instance
     * @param DcDetector DC Detector geometry
     * @return           resulting fitted cluster of hits
     */
    public List<FittedCluster> FindHitBasedClusters(List<Hit> allhits, ClusterCleanerUtilities ct,
            ClusterFitter cf, DCGeant4Factory DcDetector, boolean parallel) {

        // Create clusters and prune noise hits
        this.fillHitArray(allhits, 0);
        List<Cluster> clusters = this.findClumps(allhits, ct);
        allhits.clear();

        for (Cluster clus : clusters) {
            Collections.sort(clus);
            allhits.addAll(ct.hitListPruner(clus));
        }

        this.fillHitArray(allhits, 0);
        clusters.clear();
        clusters = this.findClumps(allhits, ct);

        // Create cluster list to be fitted
        List<FittedCluster> clusList = new ArrayList<>();
        for (Cluster c : clusters) {
            FittedCluster fc = createFittedCluster(c, ct);
            if (fc != null) clusList.add(fc);
        }

        // Fit and Split the clusters
        List<FittedCluster> fittedClusList;
        if (!parallel) fittedClusList = getFitClustersSeq(clusList, cf, ct, DcDetector);
        else           fittedClusList = getFitClustersCPUPar(clusList, DcDetector);

        // Update the clusters order and ids
        Collections.sort(fittedClusList);
        for (int ci = 0; ci < fittedClusList.size(); ++ci) fittedClusList.get(ci).set_Id(ci+1);

        return fittedClusList;
    }

    private List<FittedCluster> getFitClustersSeq(List<FittedCluster> cl, ClusterFitter cf,
            ClusterCleanerUtilities ct, DCGeant4Factory DcDetector) {

        List<FittedCluster> fcl = new ArrayList<>();
        for (FittedCluster c : cl) {
            if (fitCluster(c, cf)) fcl.addAll(ct.clusterSplitter(c, cf));
            else                   fcl.add(c);
        }

        for (FittedCluster c : fcl) {
            if (c == null) {
                fcl.remove(c);
                continue;
            }
            if (!refitCluster(c, cf, DcDetector)) fcl.remove(c);
        }

        return fcl;
    }
    private List<FittedCluster> getFitClustersCPUPar(List<FittedCluster> cl,
            DCGeant4Factory DcDetector) {

        List<FittedCluster> fcl = new CopyOnWriteArrayList<>();
        cl.parallelStream().forEach((c) -> {
            ClusterCleanerUtilities ct = new ClusterCleanerUtilities();
            ClusterFitter           cf = new ClusterFitter();
            if (fitCluster(c, cf)) fcl.addAll(ct.clusterSplitter(c, cf));
            else                   fcl.add(c);
        });

        fcl.parallelStream().forEach((c) -> {
            if (c == null) {
                fcl.remove(c);
                return;
            }

            ClusterFitter cf = new ClusterFitter();
            if (!refitCluster(c, cf, DcDetector)) fcl.remove(c);
        });

        return fcl;
    }
    private FittedCluster createFittedCluster(Cluster c, ClusterCleanerUtilities ct) {
        if (c.size() < Constants.DC_MIN_NLAYERS) return null; // TODO: Compare time with & without this line

        FittedCluster fc = new FittedCluster(c);
        ct.outOfTimersRemover(fc, true); // remove out-of-timers

        if (fc.size() < Constants.DC_MIN_NLAYERS) return null;
        return fc;
    }
    private boolean fitCluster(FittedCluster c, ClusterFitter cf) {
        cf.setFitArray(c, false);
        cf.fit(c, true);

        if      (c.get_fitProb() > Constants.HITBASEDTRKGMINFITHI2PROB) return false;
        else if (c.size() < Constants.HITBASEDTRKGNONSPLITTABLECLSSIZE) return false;
        else                                                            return true;
    }
    private boolean refitCluster(FittedCluster c, ClusterFitter cf, DCGeant4Factory DcDetector) {
        if (c == null || c.size() <= 3) return false;

        // update the hits
        for (FittedHit h : c) {
            h.set_TrkgStatus(0);
            h.updateHitPosition(DcDetector);
            h.set_AssociatedClusterID(c.get_Id());
        }

        cf.setFitArray(c, true);
        cf.fit(c, true);
        cf.setResidualDerivedParams(c, false, false, DcDetector);

        // cf.setFitArray(clus, true);
        cf.fit(c, false);
        cf.setSegmentLineParameters(c.get(0).get_Z(), c);

        if (c == null) return false;
        else           return true;
    }

    /**
     * @param fhits
     * @param tab
     * @param DcDetector DC Detector geometry
     * @param tde
     * @return
     */
    private List<FittedCluster> RecomposeClusters(List<FittedHit> fhits,
                                                  IndexedTable tab,
                                                  DCGeant4Factory DcDetector,
                                                  TimeToDistanceEstimator tde) {

        List<FittedCluster> clusters = new ArrayList<FittedCluster>();
        int NbClus = -1;
        for (FittedHit hit : fhits) {
            if (hit.get_AssociatedClusterID() == -1)    continue;
            if (hit.get_AssociatedClusterID() > NbClus) NbClus = hit.get_AssociatedClusterID();
        }

        FittedHit[][] HitArray = new FittedHit[fhits.size()][NbClus + 1];

        int index = 0;
        for (FittedHit hit : fhits) {
            if (hit.get_AssociatedClusterID() == -1) continue;
            HitArray[index][hit.get_AssociatedClusterID()] = hit;
            hit.updateHitPosition(DcDetector);

            index++;
        }

        for (int c = 0; c < NbClus + 1; c++) {
            List<FittedHit> hitlist = new ArrayList<FittedHit>();
            for (int i = 0; i < index; i++) {
                if (HitArray[i][c] != null) hitlist.add(HitArray[i][c]);
            }
            if (hitlist.size() > 0) {
                Cluster cluster = new Cluster(hitlist.get(0).get_Sector(),
                                              hitlist.get(0).get_Superlayer(),
                                              c);
                FittedCluster fcluster = new FittedCluster(cluster);
                fcluster.addAll(hitlist);
                clusters.add(fcluster);
            }
        }

        for (FittedCluster clus : clusters) {
            if (clus != null) {
                // update the hits
                for (FittedHit fhit : clus) {
                    fhit.set_TrkgStatus(0);
                    fhit.updateHitPositionWithTime(1, fhit.getB(), tab, DcDetector, tde);
                    fhit.set_AssociatedClusterID(clus.get_Id());
                    fhit.set_AssociatedHBTrackID(clus.get(0).get_AssociatedHBTrackID());
                }
            }
        }

        return clusters;
    }

    public List<FittedCluster> FindTimeBasedClusters(List<FittedHit> fhits,
                                                     ClusterFitter cf,
                                                     ClusterCleanerUtilities ct,
                                                     IndexedTable tab,
                                                     DCGeant4Factory DcDetector,
                                                     TimeToDistanceEstimator tde) {

        List<FittedCluster> clusters = new ArrayList<FittedCluster>();
        List<FittedCluster> rclusters = RecomposeClusters(fhits, tab, DcDetector, tde);

        for (FittedCluster clus : rclusters) {
            // Clean the the clusters
            FittedCluster cleanClus = ct.secondariesRemover(clus, cf, tab, DcDetector, tde);
            clus = cleanClus;

            if (clus == null) continue;

            FittedCluster LRresolvClus = ct.LRAmbiguityResolver(clus, cf, tab, DcDetector, tde);
            clus = LRresolvClus;
            if (clus == null) continue;

            // Resolve segments where there are only single hits in layers thereby resulting in a
            //     two-fold LR ambiguity hence there are 2 solutions to the segments
            int[] SumLn = new int[6];
            for (FittedHit fhit : clus) SumLn[fhit.get_Layer() - 1]++;

            boolean tryOtherClus = true;
            for (int l = 0; l < 6; l++) {
                if (SumLn[l] > 1) {
                    tryOtherClus = false;
                    break;
                }
            }

            if (tryOtherClus) {
                FittedCluster clus2 = new FittedCluster(clus.getBaseCluster());
                for (FittedHit hit : clus) {
                    if (hit.get_LeftRightAmb() == 0) continue;

                    FittedHit newhit = new FittedHit(hit.get_Sector(), hit.get_Superlayer(),
                                                     hit.get_Layer(),  hit.get_Wire(),
                                                     hit.get_TDC(), hit.get_Id());
                    newhit.set_Doca(hit.get_Doca());
                    newhit.set_DocaErr(hit.get_DocaErr());
                    newhit.setT0(hit.getT0());
                    newhit.set_Beta(hit.get_Beta());
                    newhit.setB(hit.getB());
                    newhit.set_DeltaTimeBeta(hit.get_DeltaTimeBeta());
                    newhit.setTStart(hit.getTStart());
                    newhit.setTProp(hit.getTProp());
                    newhit.setTFlight(hit.getTFlight());
                    newhit.set_Time(hit.get_Time());
                    newhit.set_Id(hit.get_Id());
                    newhit.set_TrkgStatus(hit.get_TrkgStatus());
                    newhit.set_LeftRightAmb(-hit.get_LeftRightAmb());
                    newhit.calc_CellSize(DcDetector);

                    // Assume the track angle is parallel to the layer
                    newhit.updateHitPositionWithTime(1, hit.getB(), tab, DcDetector, tde);
                    newhit.set_AssociatedClusterID(hit.get_AssociatedClusterID());
                    newhit.set_AssociatedHBTrackID(hit.get_AssociatedHBTrackID());

                    clus2.add(newhit);
                }
                cf.setFitArray(clus2, true);
                cf.fit(clus2, true);

                if (Math.abs(clus.get_Chisq() - clus2.get_Chisq()) < 1) clusters.add(clus2);
            }
            clusters.add(clus);
        }

        for (FittedCluster clus : clusters) {
            cf.setFitArray(clus, true);
            cf.fit(clus, true);

            double cosTrkAngle = 1. / Math.sqrt(1. + clus.get_clusterLineFitSlope() *
                                                     clus.get_clusterLineFitSlope());

            // Update the hits
            for (FittedHit fhit : clus) {
                fhit.updateHitPositionWithTime(cosTrkAngle, fhit.getB(), tab, DcDetector, tde);
            }

            // Iterate til convergence of trkAngle
            double Chi2Diff = 1;
            double prevChi2 = 999999999;
            double cosTrkAngleFinal = 0;
            while (Chi2Diff > 0) {
                cf.setFitArray(clus, true);
                cf.fit(clus, true);
                Chi2Diff = prevChi2 - clus.get_Chisq();
                if (Chi2Diff > 0) {
                    cosTrkAngle = 1. / Math.sqrt(1. + clus.get_clusterLineFitSlope() *
                                                      clus.get_clusterLineFitSlope());
                    // Update the hits
                    for (FittedHit fhit : clus) {
                        fhit.updateHitPositionWithTime(cosTrkAngle, fhit.getB(), tab, DcDetector, tde);
                    }
                    cosTrkAngleFinal = cosTrkAngle;
                }
                prevChi2 = clus.get_Chisq();
            }

            // Update to MP
            // calcTimeResidual = false, resetLRAmbig = false
            cf.setResidualDerivedParams(clus, false, false, DcDetector);

            for (FittedHit fhit : clus) {
                fhit.updateHitPositionWithTime(cosTrkAngleFinal, fhit.getB(), tab, DcDetector, tde);
            }

            cf.setFitArray(clus, true);
            cf.fit(clus, true);
            // calcTimeResidual = false, resetLRAmbig = false
            cf.setResidualDerivedParams(clus, true, false, DcDetector);

            cf.setFitArray(clus, true);
            cf.fit(clus, false);

            cf.setSegmentLineParameters(clus.get(0).get_Z(), clus);
        }

        return clusters;
    }

    /**
     * Selects cut to pass a hit (Not yet implemented, for now passes all hits).
     * @param hit the hit
     * @return    a flag determining if the hit is passed or not
     */
    public boolean passHitSelection(Hit hit) {
        return true;
    }

    /**
     * @param fclusters
     * @param allhits
     * @param ct
     * @param cf
     * @param event     Evio data event
     * @return
     */
    public EvioDataBank getLayerEfficiencies(List<FittedCluster> fclusters,
                                             List<Hit> allhits,
                                             ClusterCleanerUtilities ct,
                                             ClusterFitter cf,
                                             EvioDataEvent event) {

        ArrayList<Hit> clusteredHits = new ArrayList<Hit>();
        for (FittedCluster fclus : fclusters) {
            for (int k = 0; k < fclus.size(); k++) clusteredHits.add(fclus.get(k));
        }

        int[][][] EffArray = new int[6][6][6]; // 6 sectors, 6 superlayers, 6 layers
        for (int i = 0; i < 6; i++) {
            for (int j = 0; j < 6; j++) {
                for (int k = 0; k < 6; k++) {
                    EffArray[i][j][k] = -1;
                }
            }
        }

        for (int rejLy = 1; rejLy <= 6; rejLy++) {
            // Fill array of hit
            this.fillHitArray(clusteredHits, rejLy);

            // Find clumps of hits
            List<Cluster> clusters = this.findClumps(clusteredHits, ct);

            // Create cluster list to be fitted
            List<FittedCluster> selectedClusList = new ArrayList<FittedCluster>();

            for (Cluster clus : clusters) {
                FittedCluster fclus = new FittedCluster(clus);
                selectedClusList.add(fclus);
            }

            for (FittedCluster clus : selectedClusList) {
                if (clus == null) continue;
                int status = 0;

                // fit
                cf.setFitArray(clus, false);
                cf.fit(clus, true);

                for (Hit hit : allhits) {

                    if (hit.get_Sector() != clus.get_Sector()
                            || hit.get_Superlayer() != clus.get_Superlayer()
                            || hit.get_Layer() != rejLy) {
                        continue;
                    }

                    double locX = hit.calcLocY(hit.get_Layer(), hit.get_Wire());
                    double locZ = hit.get_Layer();

                    double calc_doca = Math.abs(locX - clus.get_clusterLineFitSlope() * locZ
                                                - clus.get_clusterLineFitIntercept());

                    // 2 * Math.tan(Math.PI / 6.) = 1.1547005383792515
                    if (calc_doca < 1.1547005383792515) {
                        // found a hit close enough to the track to assume that the layer is live
                        status = 1;
                    }
                    int sec = clus.get_Sector() - 1;
                    int slay = clus.get_Superlayer() - 1;
                    int lay = rejLy - 1;

                    EffArray[sec][slay][lay] = status;
                }
            }
        }

        // Now fill the bank
        int bankSize = 6 * 6 * 6;
        EvioDataBank bank = (EvioDataBank) event.getDictionary().createBank("HitBasedTrkg::LayerEffs",
                                                                            bankSize);
        int bankEntry = 0;
        for (int i = 0; i < 6; i++) {
            for (int j = 0; j < 6; j++) {
                for (int k = 0; k < 6; k++) {
                    bank.setInt("sector",     bankEntry, i + 1);
                    bank.setInt("superlayer", bankEntry, j + 1);
                    bank.setInt("layer",      bankEntry, k + 1);
                    bank.setInt("status",     bankEntry, EffArray[i][j][k]);
                    bankEntry++;
                }
            }
        }

        return bank;
    }
}
