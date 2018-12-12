package org.jlab.rec.dc.banks;

import java.util.ArrayList;
import java.util.List;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;

import org.jlab.geom.prim.Point3D;
import org.jlab.geom.prim.Vector3D;
import org.jlab.rec.dc.hit.FittedHit;
import org.jlab.rec.dc.cluster.FittedCluster;
import org.jlab.rec.dc.segment.Segment;
import org.jlab.rec.dc.cross.Cross;
import org.jlab.rec.dc.track.Track;

/**
 * A class to retrieve data from the reconstructed DC banks.
 * @author benkel
 */
public class RecoBankReader {

    // TODO: Write variation methods to be called by the DCKF engine that only
    //       return the required data to run the CrossListFinder and the
    //       TrackCandListFinder methods to avoid pulling more data than what's
    //       explicitly required.

    boolean debug = false;

    /**
     * Gets one fitted hit from a hits databank given by an index.
     * @param hBank hits bank
     * @param idx   hit's address
     * @return      the fitted hit retrieved from the bank
     */
    public FittedHit getHit(DataBank bank, int idx) {
        if (bank == null) {
            if (debug) System.out.println("[RecoBankReader.getHit] ERROR: Hits bank is null.");
            return null;
        }
        if (bank.rows() == 0) {
            if (debug) System.out.println("[RecoBankReader.getHit] ERROR: Hits bank is empty.");
            return null;
        }
        if (idx >= bank.rows()) {
            if (debug) System.out.println("[RecoBankReader.getHit] ERROR: ID given (" + idx
                                          + ") is greater than hits bank's size (" + bank.rows()
                                          + ").");
            return null;
        }

        // Read the hit
        FittedHit hit = new FittedHit((int) bank.getByte ("sector",     idx),
                                      (int) bank.getByte ("superlayer", idx),
                                      (int) bank.getByte ("layer",      idx),
                                      (int) bank.getShort("wire",       idx),
                                            bank.getInt  ("TDC",        idx),
                                      (int) bank.getShort("id",         idx));

        hit.set_Id                 ((int)    bank.getShort("id",        idx));
        hit.set_TrkgStatus         ((int)    bank.getShort("status",    idx));
        hit.set_DocaErr            ((double) bank.getFloat("docaError", idx));
        hit.set_ClusFitDoca        ((double) bank.getFloat("trkDoca",   idx));
        hit.set_lX                 ((double) bank.getFloat("LocX",      idx));
        hit.set_lY                 ((double) bank.getFloat("LocY",      idx));
        hit.set_X                  ((double) bank.getFloat("X",         idx));
        hit.set_Z                  ((double) bank.getFloat("Z",         idx));
        hit.set_LeftRightAmb       ((int)    bank.getByte ("LR",        idx));
        hit.set_AssociatedClusterID((int)    bank.getShort("clusterID", idx));
        hit.set_AssociatedHBTrackID((int)    bank.getByte ("trkID",     idx));
        hit.setB                   ((double) bank.getFloat("B",         idx));
        hit.setTProp               ((double) bank.getFloat("TProp",     idx));
        hit.setTFlight             ((double) bank.getFloat("TFlight",   idx));

        return hit;
    }

    /**
     * Gets a fitted cluster from a clusters databank given by an index, along with its referenced
     * list of fitted hits.
     * @param clBank clusters bank
     * @param hits   list of fitted hits
     * @param idx    address of the cluster
     * @return       the cluster retrieved from the bank
     */
    public FittedCluster getCluster(DataBank clBank, List<FittedHit> hits, int idx) {
        if (clBank == null) {
            if (debug) System.out.println("[RecoBankReader.getCluster] ERROR: Clusters bank is "
                                          + "null.");
            return null;
        }
        if (clBank.rows() == 0) {
            if (debug) System.out.println("[RecoBankReader.getCluster] ERROR: Clusters bank is "
                                          + "empty.");
            return null;
        }
        if (idx >= clBank.rows()) {
            if (debug) System.out.println("[RecoBankReader.getCluster] ERROR: Index ("
                                          + idx + ") is greater than clusters bank's size ("
                                          + clBank.rows() + ").");
            return null;
        }
        if (hits == null) {
            if (debug) System.out.println("[RecoBankReader.getCluster] ERROR: Hits list is "
                                          + "null.");
            return null;
        }
        if (hits.size() == 0) {
            if (debug) System.out.println("[RecoBankReader.getCluster] ERROR: Hits list is "
                                          + "empty.");
            return null;
        }

        ArrayList<Integer> hIDs    = new ArrayList<Integer>();
        ArrayList<FittedHit> fHits = new ArrayList<FittedHit>();

        for (int i = 1; i < clBank.getByte("size", idx) + 1 && i < 13; i++) {
            hIDs.add((int) clBank.getShort("Hit" + i + "_ID", idx));
        }

        // Find the list of hits from which the cluster is formed
        for (int hID : hIDs) {
            FittedHit fHit = null;
            for (FittedHit hit : hits) {
                if (hID == hit.get_Id()) {
                    fHit = hit;
                    break;
                }
            }
            if (fHit == null) {
                if (debug) System.out.println("[RecoBankReader.getCluster] ERROR: Hit with ID "
                                              + hID + " was not found in list of hits.");
                return null;
            }
            fHits.add(fHit);
        }

        FittedCluster cluster = new FittedCluster((int) clBank.getByte ("sector",       idx),
                                                  (int) clBank.getByte ("superlayer",   idx),
                                                  (int) clBank.getShort("id",           idx),
                                                  fHits);

        cluster.set_clusterLineFitSlope       ((double) clBank.getFloat("fitSlope",     idx));
        cluster.set_clusterLineFitSlopeErr    ((double) clBank.getFloat("fitSlopeErr",  idx));
        cluster.set_clusterLineFitIntercept   ((double) clBank.getFloat("fitInterc",    idx));
        cluster.set_clusterLineFitInterceptErr((double) clBank.getFloat("fitIntercErr", idx));
        cluster.set_fitProb                   ((double) clBank.getFloat("fitChisqProb", idx));
        // cluster.set_Chisq                     ((double) clBank.getFloat("fitChisqProb", idx));

        /* WARNING:
        CLUSTER DATA IN BANK CURRENTLY UNUSED:
        bank.setShort("status", i, (short) status);
        bank.setFloat("avgWire", i, (float) cluslist.get(i).getAvgwire());
        */

        return cluster;
    }

    /**
     * Gets a segment from a segments databank given by an index along with its referenced cluster.
     * @param sBank    segments bank
     * @param clusters list of clusters
     * @param idx      address of the segment
     * @return         the segment retrieved from the bank
     */
    public Segment getSegment(DataBank sBank, List<FittedCluster> clusters, int idx) {
        if (sBank == null) {
            if (debug) System.out.println("[RecoBankReader.getSegment] ERROR: Segments bank is "
                                          + "null.");
            return null;
        }
        if (sBank.rows() == 0) {
            if (debug) System.out.println("[RecoBankReader.getSegment] ERROR: Segments bank is "
                                          + "empty.");
            return null;
        }
        if (idx >= sBank.rows()) {
            if (debug) System.out.println("[RecoBankReader.getSegment] ERROR: Index ("
                                          + idx + ") is greater than segments bank's size ("
                                          + sBank.rows() + ").");
            return null;
        }
        if (clusters == null) {
            if (debug) System.out.println("[RecoBankReader.getSegment] ERROR: Cluster list "
                                          + "is null.");
            return null;
        }
        if (clusters.size() == 0 || clusters == null) {
            if (debug) System.out.println("[RecoBankReader.getSegment] ERROR: Cluster list is"
                                          + " empty.");
            return null;
        }

        // Find the cluster from which this segment is formed
        int clID = (int) sBank.getShort("Cluster_ID", idx);
        FittedCluster fCluster = null;

        for (FittedCluster cluster : clusters) {
            if (clID == cluster.get_Id()) {
                fCluster = cluster;
                break;
            }
        }
        if (fCluster == null) {
            if (debug) System.out.println("[RecoBankReader.getSegment] ERROR: Cluster with ID "
                                          + clID + " was not found in list of clusters.");
            return null;
        }

        Segment segment = new Segment(fCluster);

        double[] endPoints = {(double) sBank.getFloat("SegEndPoint1X", idx),
                              (double) sBank.getFloat("SegEndPoint1Z", idx),
                              (double) sBank.getFloat("SegEndPoint2X", idx),
                              (double) sBank.getFloat("SegEndPoint2Z", idx)};

        segment.set_SegmentEndPoints(endPoints);

        segment.get_fittedCluster()
               .set_clusterLineFitSlope((double) sBank.getFloat("fitSlope", idx));
        segment.get_fittedCluster()
               .set_clusterLineFitSlopeErr((double) sBank.getFloat("fitSlopeErr", idx));
        segment.get_fittedCluster()
               .set_clusterLineFitInterceptErr((double) sBank.getFloat("fitInterc", idx));
        segment.get_fittedCluster()
               .set_clusterLineFitInterceptErr((double) sBank.getFloat("fitIntercErr", idx));

        return segment;
    }

    /**
     * Gets a cross from a crosses databank given by an index along with its two referenced segments.
     * @param crBank   crosses bank
     * @param segments list of segments
     * @param idx      index of the cross to be retrieved
     * @return         the retrieved cross
     */
    public Cross getCross(DataBank crBank, List<Segment> segments, int idx) {
        if (crBank == null) {
            if (debug) System.out.println("[RecoBankReader.getCross] ERROR: Crosses bank is null.");
            return null;
        }
        if (crBank.rows() == 0) {
            if (debug) System.out.println("[RecoBankReader.getCross] ERROR: Crosses bank is empty.");
            return null;
        }
        if (idx >= crBank.rows()) {
            if (debug) System.out.println("[RecoBankReader.getCross] ERROR: Index given ("
                                          + idx + ") is greater than crosses bank's size ("
                                          + crBank.rows() + ").");
            return null;
        }
        if (segments == null) {
            if (debug) System.out.println("[RecoBankReader.getCross] ERROR: Segment list is "
                                          + "null.");
            return null;
        }
        if (segments.size() == 0) {
            if (debug) System.out.println("[RecoBankReader.getCross] ERROR: Segment list is"
                                          + " empty.");
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
        for (int i = 1; i < 3; i++) {
            int sID = (int) crBank.getShort("Segment" + i + "_ID", idx);
            Segment fSegment = null;
            for (Segment segment : segments) {
                if (sID == segment.get_Id()) {
                    fSegment = segment;
                    break;
                }
            }
            if (fSegment == null) {
                if (debug) System.out.println("[RecoBankReader.getCross] ERROR: Segment with ID "
                                              + sID + " was not found in list of segments.");
                return null;
            }

            if (i == 1)      cross.set_Segment1(fSegment);
            else if (i == 2) cross.set_Segment2(fSegment);
            cross.add(fSegment);
        }

        if (cross.get_Segment1().get_Id() == -1 || cross.get_Segment2().get_Id() == -1)
            cross.isPseudoCross = true;

        cross.set_CrossDirIntersSegWires();

        return cross;
    }

    /**
     * Gets a track from a databank, given by an index, along with its three referenced segments.
     * @param trBank  tracks bank
     * @param crosses list of crosses
     * @param idx     index of the track
     * @return        the retrieved track
     */
    public Track getHBTrack(DataBank trBank, List<Cross> crosses, int idx) {
        if (trBank == null) {
            if (debug) System.out.println("[RecoBankReader.getHBTrack] ERROR: Tracks bank is "
                                          + "null.");
            return null;
        }
        if (trBank.rows() == 0) {
            if (debug) System.out.println("[RecoBankReader.getHBTrack] ERROR: Tracks bank is "
                                          + "empty.");
            return null;
        }
        if (idx > trBank.rows()) {
            if (debug) System.out.println("[RecoBankReader.getHBTrack] ERROR: Index given (" +
                                          + idx + ") is greater than tracks bank's size (" +
                                          trBank.rows() + ").");
            return null;
        }
        if (crosses == null) {
            if (debug) System.out.println("[RecoBankReader.getHBTrack] ERROR: Cross list given is "
                                          + "null.");
            return null;
        }
        if (crosses.size() == 0) {
            if (debug) System.out.println("[RecoBankReader.getHBTrack] ERROR: Cross list given is "
                                          + "empty.");
            return null;
        }

        Track track = new Track();
        track.set_Id    ((int) trBank.getShort("id",     idx));
        track.set_Sector((int) trBank.getByte ("sector", idx));
        track.set_Q     ((int) trBank.getByte ("q",      idx));

        // TODO: The status is changed when being written into the bank, see if this affects anything.
        track.set_Status((int) trBank.getShort("status", idx));

        for (int i = 1; i < 4; i++) {
            int cID = (int) trBank.getShort("Cross" + i + "_ID", idx);
            Cross fCross = null;
            for (Cross cross : crosses) {
                if (cID == cross.get_Id()) {
                    fCross = cross;
                    break;
                }
            }
            if (fCross == null) {
                if (debug) System.out.println("[RecoBankReader.getHBTrack] ERROR: Cross with ID "
                                              + cID + " was not found in cross list.");
                return null;
            }

            track.add(fCross);
        }

        if (trBank.getFloat("c1_x") != null) {
            track.set_PreRegion1CrossPoint(new Point3D((double) trBank.getFloat("c1_x", idx),
                                                       (double) trBank.getFloat("c1_y", idx),
                                                       (double) trBank.getFloat("c1_z", idx)));
            track.set_PreRegion1CrossDir  (new Point3D((double) trBank.getFloat("c1_ux", idx),
                                                       (double) trBank.getFloat("c1_uy", idx),
                                                       (double) trBank.getFloat("c1_uz", idx)));
        }
        if (trBank.getFloat("c3_x") != null) {
            track.set_PostRegion3CrossPoint(new Point3D((double) trBank.getFloat("c3_x", idx),
                                                        (double) trBank.getFloat("c3_y", idx),
                                                        (double) trBank.getFloat("c3_z", idx)));
            track.set_PostRegion3CrossDir  (new Point3D((double) trBank.getFloat("c3_ux", idx),
                                                        (double) trBank.getFloat("c3_uy", idx),
                                                        (double) trBank.getFloat("c3_uz", idx)));
        }
        if (trBank.getFloat("t1_x") != null) {
            track.set_Region1TrackX(new Point3D((double) trBank.getFloat("t1_x", idx),
                                                (double) trBank.getFloat("t1_y", idx),
                                                (double) trBank.getFloat("t1_z", idx)));
            track.set_Region1TrackP(new Point3D((double) trBank.getFloat("t1_px", idx),
                                                (double) trBank.getFloat("t1_py", idx),
                                                (double) trBank.getFloat("t1_pz", idx)));
        }

        track.set_TotPathLen          ((double) trBank.getFloat("pathlength", idx));
        track.set_Vtx0(new Point3D    ((double) trBank.getFloat("Vtx0_x", idx),
                                       (double) trBank.getFloat("Vtx0_y", idx),
                                       (double) trBank.getFloat("Vtx0_z", idx)));
        track.set_pAtOrig(new Vector3D((double) trBank.getFloat("p0_x", idx),
                                       (double) trBank.getFloat("p0_y", idx),
                                       (double) trBank.getFloat("p0_z", idx)));
        track.set_FitChi2             ((double) trBank.getFloat("chi2", idx));
        track.set_FitNDF              ((int)    trBank.getShort("ndf", idx));

        return track;
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
    public static void printFullInfo(List<Cross> cList) {
        System.out.println("Crosses info:");
        for (int c = 0; c < cList.size(); c++) {
            System.out.println(cList.get(c).getDetailedInfo());

            System.out.println(cList.get(c).get_Segment1().getDetailedInfo());
            System.out.println(cList.get(c).get_Segment2().getDetailedInfo());

            System.out.println(cList.get(c).get_Segment1().get_fittedCluster().getDetailedInfo());
            System.out.println(cList.get(c).get_Segment2().get_fittedCluster().getDetailedInfo());

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
