package org.jlab.rec.dc.banks;

import java.util.ArrayList;
import java.util.List;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;

import org.jlab.geom.prim.Point3D;
import org.jlab.rec.dc.hit.FittedHit;
import org.jlab.rec.dc.cluster.FittedCluster;
import org.jlab.rec.dc.segment.Segment;
import org.jlab.rec.dc.cross.Cross;
import org.jlab.rec.dc.track.Track;

/**
 * A class to retrieve the data from the reconstructed DC banks.
 * @author benkel
 */
public class RecoBankReader {

    // TODO: Write a manager method to avoid creating multiple copies of
    //       referenced objects.
    // TODO: Write variation methods to be called by the DCKF engine that only
    //       return the required data to run the CrossListFinder and the
    //       TrackCandListFinder methods to avoid pulling more data than what's
    //       explicitly required.

    /**
     * Gets one fitted hit from a hits databank given by an index.
     * @param hBank bank containing the fitted hits
     * @param addr  address of the hit to be retrieved
     * @return      the hit retrieved from the bank
     */
    public FittedHit getHit(DataBank bank, int addr) {
        if (bank.rows() == 0) {
            System.out.println("[DCKF] ERROR: hits bank is empty.");
            return null;
        }
        if (bank.rows() == 0 || addr >= bank.rows()) {
            System.out.println("[DCKF] ERROR: id given is greater than hits " +
                               "bank's size.");
            return null;
        }

        FittedHit hit = new FittedHit((int) bank.getByte ("sector",     addr),
                                      (int) bank.getByte ("superlayer", addr),
                                      (int) bank.getByte ("layer",      addr),
                                      (int) bank.getShort("wire",       addr),
                                            bank.getInt  ("TDC",        addr),
                                      (int) bank.getShort("id",         addr));

        hit.set_Id                 ((int)    bank.getShort("id",        addr));
        hit.set_TrkgStatus         ((int)    bank.getShort("status",    addr));
        hit.set_DocaErr            ((double) bank.getFloat("docaError", addr));
        hit.set_ClusFitDoca        ((double) bank.getFloat("trkDoca",   addr));
        hit.set_lX                 ((double) bank.getFloat("LocX",      addr));
        hit.set_lY                 ((double) bank.getFloat("LocY",      addr));
        hit.set_X                  ((double) bank.getFloat("X",         addr));
        hit.set_Z                  ((double) bank.getFloat("Z",         addr));
        hit.set_LeftRightAmb       ((int)    bank.getByte ("LR",        addr));
        hit.set_AssociatedClusterID((int)    bank.getShort("clusterID", addr));
        hit.set_AssociatedHBTrackID((int)    bank.getByte ("trkID",     addr));
        hit.setB                   ((double) bank.getFloat("B",         addr));
        hit.setTProp               ((double) bank.getFloat("TProp",     addr));
        hit.setTFlight             ((double) bank.getFloat("TFlight",   addr));

        return hit;
    }

    /**
     * Gets one fitted cluster from a clusters databank given by an index along with its referenced
     * list of hits.
     * @param clBank bank containing the fitted clusters
     * @param hBank  bank containing the fitted hits, referenced by the fitted clusters
     * @param clAddr address of the cluster to be retrieved
     * @return       the cluster retrieved from the bank
     */
    public FittedCluster getCluster(DataBank clBank, DataBank hBank, int clAddr) {
        if (clBank.rows() == 0) {
            System.out.println("[DCKF] ERROR: clusters bank is empty.");
            return null;
        }
        if (clAddr >= clBank.rows()) {
            System.out.println("[DCKF] ERROR: index given is great than clusters bank's size.");
            return null;
        }

        ArrayList<Integer> hIDs    = new ArrayList<Integer>();
        ArrayList<FittedHit> fHits = new ArrayList<FittedHit>();

        for (int i = 1; i < clBank.getByte("size", clAddr) + 1; i++) {
            hIDs.add((int) clBank.getShort("Hit" + i + "_ID", clAddr));
        }

        int nHits = hBank.rows();
        if (nHits == 0) {
            System.out.println("[DCKF] ERROR: hits bank is empty.");
            return null;
        }

        for (int hID : hIDs) {
            int hAddress = -1;
            for (int i = 0; i < nHits; i++) {
                if (hID == (int) hBank.getShort("id", i)) {
                    hAddress = i;
                    break;
                }
            }
            if (hAddress == -1) {
                System.out.println("[DCKF] ERROR: hit id was not found in " +
                                   "hits bank");
                return null;
            }
            FittedHit fHit = getHit(hBank, hAddress);
            fHits.add(fHit);
        }

        FittedCluster cluster = new FittedCluster((int) clBank.getByte ("sector",       clAddr),
                                                  (int) clBank.getByte ("superlayer",   clAddr),
                                                  (int) clBank.getShort("id",           clAddr),
                                                  fHits);

        cluster.set_clusterLineFitSlope       ((double) clBank.getFloat("fitSlope",     clAddr));
        cluster.set_clusterLineFitSlopeErr    ((double) clBank.getFloat("fitSlopeErr",  clAddr));
        cluster.set_clusterLineFitIntercept   ((double) clBank.getFloat("fitInterc",    clAddr));
        cluster.set_clusterLineFitInterceptErr((double) clBank.getFloat("fitIntercErr", clAddr));
        cluster.set_fitProb                   ((double) clBank.getFloat("fitChisqProb", clAddr));
        // cluster.set_Chisq                     ((double) clBank.getFloat("fitChisqProb", clAddr));



        /* WARNING:
        CLUSTER DATA IN BANK CURRENTLY UNUSED:
        bank.setShort("status", i, (short) status);
        bank.setFloat("avgWire", i, (float) cluslist.get(i).getAvgwire());
        bank.setFloat("fitChisqProb", i, (float) ProbChi2perNDF.prob(chi2, cluslist.get(i).size() - 2));
        */

        return cluster;
    }

    /**
     * Gets one segment from a segments databank given by an index along with its referenced fitted
     * cluster.
     * @param sBank  bank containing the segments
     * @param clBank bank containing the fitted clusters, referenced by the segments
     * @param hBank  bank containing the fitted hits, referenced by the fitted clusters
     * @param sAddr  address of the segment to be retrieved
     * @return       the segment retrieved from the bank
     */
    public Segment getSegment(DataBank sBank, DataBank clBank,
                              DataBank hBank, int sAddr) {
        if (sBank.rows() == 0) {
            System.out.println("[DCKF] ERROR: segments bank is empty.");
            return null;
        }
        if (sAddr >= sBank.rows()) {
            System.out.println("[DCKF] ERROR: index given is greater than " +
                               "segments bank's size.");
            return null;
        }
        int nClusters = clBank.rows();
        if (nClusters == 0) {
            System.out.println("[DCKF] ERROR: clusters bank is empty.");
            return null;
        }
        int clID = (int) sBank.getShort("Cluster_ID", sAddr);
        int clAddress = -1;
        for (int i = 0; i < nClusters; i++) {
            if (clID == (int) clBank.getShort("id", i)) {
                clAddress = i;
                break;
            }
        }
        if (clAddress == -1) {
            System.out.println("[DCKF] ERROR: cluster id was not found in " +
                               " clusters bank.");
            return null;
        }

        FittedCluster fCluster = getCluster(clBank, hBank, clAddress);

        Segment segment = new Segment(fCluster);

        double[] endPoints = {(double) sBank.getFloat("SegEndPoint1X", sAddr),
                              (double) sBank.getFloat("SegEndPoint1Z", sAddr),
                              (double) sBank.getFloat("SegEndPoint2X", sAddr),
                              (double) sBank.getFloat("SegEndPoint2Z", sAddr)};

        segment.set_SegmentEndPoints(endPoints);

        segment.get_fittedCluster()
               .set_clusterLineFitSlope((double) sBank.getFloat("fitSlope", sAddr));
        segment.get_fittedCluster()
               .set_clusterLineFitSlopeErr((double) sBank.getFloat("fitSlopeErr", sAddr));
        segment.get_fittedCluster()
               .set_clusterLineFitInterceptErr((double) sBank.getFloat("fitInterc", sAddr));
        segment.get_fittedCluster()
               .set_clusterLineFitInterceptErr((double) sBank.getFloat("fitIntercErr", sAddr));

        /* WARNING:
        SEGMENT DATA IN BANK CURRENTLY UNUSED:
        bank.setFloat("avgWire", i, (float) cls.getAvgwire());
        bank.setByte("size", i, (byte) seglist.get(i).size());
        */
        return segment;
    }

    /**
     * Gets one cross from a crosses databank given by an index along with its two referenced
     * segments.
     * @param crBank bank containing the crosses
     * @param sBank  bank containing the segments, referenced by the crosses
     * @param clBank bank containing the fitted clusters, referenced by the segments
     * @param hBank  bank containing the fitted hits, referenced by the fitted clusters
     * @param idx    index of the cross to be retrieved
     * @return       the retrieved cross
     */
    public Cross getCross(DataBank crBank, DataBank sBank,
                          DataBank clBank, DataBank hBank,
                          int idx) {
        if (crBank.rows() == 0) {
            System.out.println("[DCKF] ERROR: crosses bank is empty.");
            return null;
        }
        if (idx > crBank.rows()) {
            System.out.println("[DCKF] ERROR: index given is greater than " +
                               "crosses bank's size.");
            return null;
        }
        Cross cross = new Cross((int)  crBank.getByte ("sector", idx),
                                (int)  crBank.getByte ("region", idx),
                                (int)  crBank.getShort("id",     idx));

        cross.set_Point   (new Point3D(crBank.getFloat("x",      idx),
                                       crBank.getFloat("y",      idx),
                                       crBank.getFloat("z",      idx)));
        cross.set_PointErr(new Point3D(crBank.getFloat("err_x",  idx),
                                       crBank.getFloat("err_y",  idx),
                                       crBank.getFloat("err_z",  idx)));
        cross.set_Dir     (new Point3D(crBank.getFloat("ux",     idx),
                                       crBank.getFloat("uy",     idx),
                                       crBank.getFloat("uz",     idx)));
        cross.set_DirErr  (new Point3D(crBank.getFloat("err_ux", idx),
                                       crBank.getFloat("err_uy", idx),
                                       crBank.getFloat("err_uz", idx)));

        // get the segments
        int nSegments = sBank.rows();
        if (nSegments == 0) {
            System.out.println("[DCKF] ERROR: segments bank is empty.");
            return null;
        }

        for (int i = 1; i < 3; i++) {
            int sID = (int) crBank.getShort("Segment" + i + "_ID", idx);
            int sAddress = -1;
            for (int j = 0; j < nSegments; j++) {
                if (sID == (int) sBank.getShort("id", j)) {
                    sAddress = j;
                    break;
                }
            }
            if (sAddress == -1) {
                System.out.println("[DCKF] ERROR: segment id was not found " +
                                   "in bank.");
                return null;
            }

            Segment segment = getSegment(sBank, clBank, hBank, sAddress);

            if (i == 1)      cross.set_Segment1(segment);
            else if (i == 2) cross.set_Segment2(segment);
            cross.add(segment);
        }

        if (cross.get_Segment1().get_Id() == -1 ||
            cross.get_Segment2().get_Id() == -1) {
            cross.isPseudoCross = true;
        }

        cross.set_CrossDirIntersSegWires();

        return cross;
    }

    /**
     * Prints the data associated to one cross along with its first segment, cluster and hit.
     * @param cross the cross to be printed
     */
    public static void printSample(Cross cross) {
        System.out.println(cross.getDetailedInfo());
        System.out.println(cross.get_Segment1().getDetailedInfo());
        System.out.println(cross.get_Segment1().get_fittedCluster().getDetailedInfo());
        System.out.println(cross.get_Segment1().get_fittedCluster().get(0).getDetailedInfo());
    }

    /**
     * Prints the data contained in a list of crosses along with all the referenced objects.
     * Warning: The output is usually very extensive.
     * @param cList the list of crosses to be printed
     */
    public static void printCrossesInfo(List<Cross> cList) {
        System.out.println("Crosses info:");
        for (int c = 0; c < cList.size(); c++) {
            System.out.println(cList.get(c).getDetailedInfo());

            System.out.println(cList.get(c).get_Segment1().getDetailedInfo());
            System.out.println(cList.get(c).get_Segment2().getDetailedInfo());

            System.out.println(cList.get(c).get_Segment1().get_fittedCluster().getDetailedInfo());
            System.out.println(cList.get(c).get_Segment2().get_fittedCluster().getDetailedInfo());

            /**************************************************************************************/
            FittedCluster cluster1 = cList.get(c).get_Segment1().get_fittedCluster();
            FittedCluster cluster2 = cList.get(c).get_Segment2().get_fittedCluster();
            System.out.println("--------------------------------------------");
            System.out.println("Hits in Cluster " + cluster1.get_Id() + ":");
            System.out.println("--------------------------------------------");
            for (int i = 0; i < cluster1.size(); i++) {
                System.out.println(cluster1.get(i).getDetailedInfo());
            }
            System.out.println("--------------------------------------------");
            System.out.println("Hits in Cluster " + cluster2.get_Id() + ":");
            System.out.println("--------------------------------------------");
            for (int i = 0; i < cluster2.size(); i++) {
                System.out.println(cluster2.get(i).getDetailedInfo());
            }
            System.out.println("--------------------------------------------");
            System.out.println("");
        }
        return;
    }
}
