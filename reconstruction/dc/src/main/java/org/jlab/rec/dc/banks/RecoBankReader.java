package org.jlab.rec.dc.banks;

import java.util.ArrayList;
import java.util.List;
import Jama.Matrix;
import org.jlab.geom.prim.Point3D;
import org.jlab.geom.prim.Vector3D;
import org.jlab.geom.prim.Line3D;
import org.jlab.geom.prim.Plane3D;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;

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

    // TODO: Write variation methods to be called by the DCKF engine that only return the required
    //       data to run the CrossListFinder and the TrackCandListFinder methods to avoid pulling
    //       more data than what's explicitly required.

    boolean debug = true;

    /**
     * Asserts that a bank and an index accessed from it are valid.
     * @param bank       bank to be validated
     * @param idx        index of the bank trying to be accessed
     * @param callerName name of the type of bank being asserted (HBhits, HBclusters, etc)
     * @return           false if the bank and the access is valid, true otherwise
     */
    public boolean validateBank(DataBank bank, int idx, String callerName) {
        if (bank == null) {
            if (debug) System.out.println("[RBR] ERROR: " + callerName + " bank is null.");
            return true;
        }
        if (bank.rows() == 0) {
            if (debug) System.out.println("[RBR] ERROR: " + callerName + " bank is empty.");
            return true;
        }
        if (idx >= bank.rows()) {
            if (debug) System.out.println("[RBR] ERROR: index given (" + idx + ") is larger than "
                                        + callerName + " bank's size (" + bank.rows() + ").");
            return true;
        }
        return false;
    }

    /**
     * Asserts that a list is valid.
     * @param list       list to be validated
     * @param listName   name of the list being asserted (hits, clusters, etc), used for reporting
     *                   errors.
     * @param callerName name of the type of bank being asserted (HBsegments, HBcrosses, etc), used
     *                   for reporting errors
     * @return           false if the list is valid, true otherwise
     */
    public boolean validateList(List<?> list, String listName, String callerName) {
        if (list == null) {
            if (debug) System.out.println("[RBR] ERROR: " + listName + " is null.");
            return true;
        }
        if (list.size() == 0) {
            if (debug) System.out.println("[RBR] ERROR: " + listName + " is empty.");
            return true;
        }
        return false;
    }

    /**
     * Gets one fitted hit from a hits databank given by an index.
     * @param hBank hits bank
     * @param idx   hit's address
     * @return      the fitted hit retrieved from the bank
     */
    public FittedHit getHit(DataBank bank, int idx) {
        if (validateBank(bank, idx, "HBhits")) return null;

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

        // TODO: v CHECK NEW STUFF v
        hit.set_CellSize              ((double) bank.getFloat("cellSize", idx));
        hit.set_XMP                   ((double) bank.getFloat("XMP", idx));
        hit.set_Residual              ((double) bank.getFloat("residual", idx));
        hit.set_TimeResidual          ((double) bank.getFloat("timeResidual", idx));
        hit.set_QualityFac            ((double) bank.getFloat("qualityFac", idx));
        hit.set_TrkgStatus            ((int)    bank.getByte ("trkStatus", idx));
        hit.set_ClusFitDoca           ((double) bank.getFloat("clusFitDoca", idx));
        hit.set_TrkFitDoca            ((double) bank.getFloat("trkFitDoca", idx));
        // hit.set_TimeToDistance        ((double) bank.getFloat("timeToDistance", idx));
        hit.set_Beta                  ((double) bank.getFloat("beta", idx));
        hit.set_Doca                  ((double) bank.getFloat("doca", idx));
        hit._lr =                      (int)    bank.getByte("lr", idx);
        // hit.setCrossDirIntersWire(new Point3D((double) bank.getFloat("crossDirIntersWireX", idx),
        //                                            (double) bank.getFloat("crossDirIntersWireY", idx),
        //                                            (double) bank.getFloat("crossDirIntersWireZ", idx)));
        // hit.setSignalPropagAlongWire  ((double) bank.getFloat("signalPropagAlongWire", idx));
        hit.setSignalTimeOfFlight     ();
        hit.setT0                     ((double) bank.getFloat("t0", idx));
        hit.setTStart                 ((double) bank.getFloat("tStart", idx));
        hit.set_Time                  ((double) bank.getFloat("time", idx));
        hit.set_OutOfTimeFlag         (bank.getByte ("outOfTimeFlag", idx) == 1 ? true : false);
        hit.set_WireLength            ((double) bank.getFloat("wireLength", idx));
        hit.set_WireMaxSag            ((double) bank.getFloat("wireMaxSag", idx));
        hit.set_TrkResid              ((double) bank.getFloat("trkResid", idx));
        hit.set_DeltaTimeBeta         ((double) bank.getFloat("deltaTimeBeta", idx));
        // TODO: new stuff is read up to here.

        return hit;
    }

    /**
     * Gets a fitted cluster from a clusters databank given by an index, along with its referenced
     * list of fitted hits.
     * @param bank   clusters bank
     * @param hits   list of fitted hits
     * @param idx    address of the cluster
     * @return       the cluster retrieved from the bank
     */
    public FittedCluster getCluster(DataBank bank, List<FittedHit> hits, int idx) {
        if (validateBank(bank, idx, "HBclusters"))    return null;
        if (validateList(hits, "hits", "HBclusters")) return null;

        ArrayList<Integer> hIDs    = new ArrayList<Integer>();
        ArrayList<FittedHit> fHits = new ArrayList<FittedHit>();

        for (int i = 1; i < bank.getByte("size", idx) + 1 && i < 13; i++) {
            hIDs.add((int) bank.getShort("Hit" + i + "_ID", idx));
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

        FittedCluster cluster = new FittedCluster((int) bank.getByte ("sector",       idx),
                                                  (int) bank.getByte ("superlayer",   idx),
                                                  (int) bank.getShort("id",           idx),
                                                  fHits);

        cluster.set_clusterLineFitSlope       ((double) bank.getFloat("fitSlope",     idx));
        cluster.set_clusterLineFitSlopeErr    ((double) bank.getFloat("fitSlopeErr",  idx));
        cluster.set_clusterLineFitIntercept   ((double) bank.getFloat("fitInterc",    idx));
        cluster.set_clusterLineFitInterceptErr((double) bank.getFloat("fitIntercErr", idx));
        cluster.set_clusterLineFitSlIntCov    ((double) bank.getFloat("fitSlIntCov",  idx));
        cluster.set_fitProb                   ((double) bank.getFloat("fitChisqProb", idx));
        cluster.set_Chisq                     ((double) bank.getFloat("chisqProb",    idx));

        // v TODO: NEW CODE, PROBABLY NOT AFFECTING ANYTHING. TRY COMMENTING IT ALL.
        cluster.set_clusLine(new Line3D((double) bank.getFloat("clusLine1X", idx),
                                        (double) bank.getFloat("clusLine1Y", idx),
                                        (double) bank.getFloat("clusLine1Z", idx),
                                        (double) bank.getFloat("clusLine2X", idx),
                                        (double) bank.getFloat("clusLine2Y", idx),
                                        (double) bank.getFloat("clusLine2Z", idx)));
        cluster.set_clusLineErr(new Line3D((double) bank.getFloat("clusLineErr1X", idx),
                                           (double) bank.getFloat("clusLineErr1Y", idx),
                                           (double) bank.getFloat("clusLineErr1Z", idx),
                                           (double) bank.getFloat("clusLineErr2X", idx),
                                           (double) bank.getFloat("clusLineErr2Y", idx),
                                           (double) bank.getFloat("clusLineErr2Z", idx)));
        cluster.set_clusterLineFitSlopeMP((double) bank.getFloat("clusterLineFitSlopeMP", idx));
        cluster.set_clusterLineFitSlopeErrMP((double) bank.getFloat("clusterLineFitSlopeErrMP", idx));
        cluster.set_clusterLineFitInterceptMP((double) bank.getFloat("clusterLineFitInterceptMP", idx));
        cluster.set_clusterLineFitInterceptErrMP((double) bank.getFloat("clusterLineFitInterceptErrMP", idx));
        // TODO: UP TO HERE

        /* WARNING:
        CLUSTER DATA IN BANK CURRENTLY UNUSED:
        bank.setShort("status", i, (short) status);
        bank.setFloat("avgWire", i, (float) cluslist.get(i).getAvgwire());
        */

        return cluster;
    }

    /**
     * Gets a segment from a segments databank given by an index along with its referenced cluster.
     * @param bank     segments bank
     * @param clusters list of clusters
     * @param idx      address of the segment
     * @return         the segment retrieved from the bank
     */
    public Segment getSegment(DataBank bank, List<FittedCluster> clusters, int idx) {
        if (validateBank(bank, idx, "HBsegments"))            return null;
        if (validateList(clusters, "clusters", "HBsegments")) return null;

        // Find the cluster from which this segment is formed
        int clID = (int) bank.getShort("Cluster_ID", idx);
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

        Segment segment = new Segment(fCluster); // TODO: Check what is given to the segment in this line!!

        double[] endPoints = {(double) bank.getFloat("SegEndPoint1X", idx),
                              (double) bank.getFloat("SegEndPoint1Z", idx),
                              (double) bank.getFloat("SegEndPoint2X", idx),
                              (double) bank.getFloat("SegEndPoint2Z", idx)};

        segment.set_SegmentEndPoints(endPoints);

        // v TODO: NEW CODE, TEST THOROUGHLY
        segment.isOnTrack = bank.getByte("isOnTrack", idx) == 1 ? true : false;
        segment.set_ResiSum((double) bank.getFloat("resiSum", idx));
        segment.set_TimeSum((double) bank.getFloat("timeSum", idx));
        segment.set_Status((int) bank.getByte("status", idx));
        segment.associatedCrossId = (int) bank.getShort("associatedCrossId", idx);

        segment.set_fitPlane(new Plane3D((double) bank.getFloat("fitPlane_px", idx),
                                         (double) bank.getFloat("fitPlane_py", idx),
                                         (double) bank.getFloat("fitPlane_pz", idx),
                                         (double) bank.getFloat("fitPlane_nx", idx),
                                         (double) bank.getFloat("fitPlane_ny", idx),
                                         (double) bank.getFloat("fitPlane_nz", idx)));
        // TODO: UP TO HERE

        // segment.get_fittedCluster()
        //        .set_clusterLineFitSlope((double) bank.getFloat("fitSlope", idx));
        // segment.get_fittedCluster()
        //        .set_clusterLineFitSlopeErr((double) bank.getFloat("fitSlopeErr", idx));
        // segment.get_fittedCluster()
        //        .set_clusterLineFitInterceptErr((double) bank.getFloat("fitInterc", idx));
        // segment.get_fittedCluster()
        //        .set_clusterLineFitInterceptErr((double) bank.getFloat("fitIntercErr", idx));

        return segment;
    }

    /**
     * Gets a cross from a crosses databank given by an index along with its two referenced segments.
     * @param bank     crosses bank
     * @param segments list of segments
     * @param idx      index of the cross to be retrieved
     * @return         the retrieved cross
     */
    public Cross getCross(DataBank bank, List<Segment> segments, int idx) {
        if (validateBank(bank, idx, "HBcrosses"))            return null;
        if (validateList(segments, "segments", "HBcrosses")) return null;

        Cross cross = new Cross((int)  bank.getByte ("sector", idx),
                                (int)  bank.getByte ("region", idx),
                                (int)  bank.getShort("id",     idx));

        cross.set_Point   (new Point3D(bank.getFloat("x",      idx),
                                       bank.getFloat("y",      idx),
                                       bank.getFloat("z",      idx)));
        cross.set_PointErr(new Point3D(bank.getFloat("err_x",  idx),
                                       bank.getFloat("err_y",  idx),
                                       bank.getFloat("err_z",  idx)));
        cross.set_Dir     (new Point3D(bank.getFloat("ux",     idx),
                                       bank.getFloat("uy",     idx),
                                       bank.getFloat("uz",     idx)));
        cross.set_DirErr  (new Point3D(bank.getFloat("err_ux", idx),
                                       bank.getFloat("err_uy", idx),
                                       bank.getFloat("err_uz", idx)));

        // get the segments
        for (int i = 1; i < 3; i++) {
            int sID = (int) bank.getShort("Segment" + i + "_ID", idx);
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
     * Gets a track from a databank, given by an index, along with its three referenced crosses.
     * @param bank    tracks bank
     * @param crosses list of crosses
     * @param idx     index of the track
     * @return        the retrieved track
     */
    public Track getHBTrack(DataBank bank, List<Cross> crosses, int idx) {
        if (validateBank(bank, idx, "HBtracks"))            return null;
        if (validateList(crosses, "crosses", "HBsegments")) return null;

        Track track = new Track();
        track.set_Id    ((int) bank.getShort("id",     idx));
        track.set_Sector((int) bank.getByte ("sector", idx));
        track.set_Q     ((int) bank.getByte ("q",      idx));

        // TODO: The status is changed when being written into the bank, see if this affects anything.
        track.set_Status((int) bank.getShort("status", idx));

        for (int i = 1; i < 4; i++) {
            int cID = (int) bank.getShort("Cross" + i + "_ID", idx);
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

        if (bank.getFloat("c1_x") != null) {
            track.set_PreRegion1CrossPoint(new Point3D((double) bank.getFloat("c1_x", idx),
                                                       (double) bank.getFloat("c1_y", idx),
                                                       (double) bank.getFloat("c1_z", idx)));
            track.set_PreRegion1CrossDir  (new Point3D((double) bank.getFloat("c1_ux", idx),
                                                       (double) bank.getFloat("c1_uy", idx),
                                                       (double) bank.getFloat("c1_uz", idx)));
        }
        if (bank.getFloat("c3_x") != null) {
            track.set_PostRegion3CrossPoint(new Point3D((double) bank.getFloat("c3_x", idx),
                                                        (double) bank.getFloat("c3_y", idx),
                                                        (double) bank.getFloat("c3_z", idx)));
            track.set_PostRegion3CrossDir  (new Point3D((double) bank.getFloat("c3_ux", idx),
                                                        (double) bank.getFloat("c3_uy", idx),
                                                        (double) bank.getFloat("c3_uz", idx)));
        }
        if (bank.getFloat("t1_x") != null) {
            track.set_Region1TrackX(new Point3D((double) bank.getFloat("t1_x", idx),
                                                (double) bank.getFloat("t1_y", idx),
                                                (double) bank.getFloat("t1_z", idx)));
            track.set_Region1TrackP(new Point3D((double) bank.getFloat("t1_px", idx),
                                                (double) bank.getFloat("t1_py", idx),
                                                (double) bank.getFloat("t1_pz", idx)));
        }

        track.set_TotPathLen          ((double) bank.getFloat("pathlength", idx));
        track.set_Vtx0(new Point3D    ((double) bank.getFloat("Vtx0_x", idx),
                                       (double) bank.getFloat("Vtx0_y", idx),
                                       (double) bank.getFloat("Vtx0_z", idx)));
        track.set_pAtOrig(new Vector3D((double) bank.getFloat("p0_x", idx),
                                       (double) bank.getFloat("p0_y", idx),
                                       (double) bank.getFloat("p0_z", idx)));
        track.set_FitChi2             ((double) bank.getFloat("chi2", idx));
        track.set_FitNDF              ((int)    bank.getShort("ndf", idx));

        // TODO: v CHECK THESE VARIABLES v
        track.set_IntegralBdl((double) bank.getFloat("_IntegralBdl", idx));
        track.set_PathLength((double) bank.getFloat("_pathLength", idx));
        track.set_P((double) bank.getFloat("_P", idx));
        track.set_Vtx0(new Point3D((double) bank.getFloat("_Vtx0_x", idx),
                                   (double) bank.getFloat("_Vtx0_y", idx),
                                   (double) bank.getFloat("_Vtx0_z", idx)));
        track.set_Vtx0(new Point3D((double) bank.getFloat("_pAtOrig_TiltedCS_x", idx),
                                   (double) bank.getFloat("_pAtOrig_TiltedCS_y", idx),
                                   (double) bank.getFloat("_pAtOrig_TiltedCS_z", idx)));

        track.fit_Successful = bank.getByte("fit_Successful", idx) == 1 ? true : false;
        track.set_MissingSuperlayer(bank.getShort("_missingSuperlayer", idx));
        track.set_FitConvergenceStatus(bank.getShort("_fitConvergenceStatus", idx));
        track.b[0] = bank.getFloat("b_0", idx);
        track.b[1] = bank.getFloat("b_1", idx);
        track.b[2] = bank.getFloat("b_2", idx);
        // TODO: UP UNTIL HERE

        return track;
    }

    /**
     * Gets a TB covariance matrix from a bank.
     * @param bank    covariance matrices bank
     * @param idx     index of the track
     * @return        the retrieved matrix
     */
    public Matrix getTBCovMat(DataBank bank, int idx) {
        if (validateBank(bank, idx, "TBCovMat")) return null;

        double[] tmpArr = {(double) bank.getFloat("C11", idx),
                           (double) bank.getFloat("C12", idx),
                           (double) bank.getFloat("C13", idx),
                           (double) bank.getFloat("C14", idx),
                           (double) bank.getFloat("C15", idx),
                           (double) bank.getFloat("C21", idx),
                           (double) bank.getFloat("C22", idx),
                           (double) bank.getFloat("C23", idx),
                           (double) bank.getFloat("C24", idx),
                           (double) bank.getFloat("C25", idx),
                           (double) bank.getFloat("C31", idx),
                           (double) bank.getFloat("C32", idx),
                           (double) bank.getFloat("C33", idx),
                           (double) bank.getFloat("C34", idx),
                           (double) bank.getFloat("C35", idx),
                           (double) bank.getFloat("C41", idx),
                           (double) bank.getFloat("C42", idx),
                           (double) bank.getFloat("C43", idx),
                           (double) bank.getFloat("C44", idx),
                           (double) bank.getFloat("C45", idx),
                           (double) bank.getFloat("C51", idx),
                           (double) bank.getFloat("C52", idx),
                           (double) bank.getFloat("C53", idx),
                           (double) bank.getFloat("C54", idx),
                           (double) bank.getFloat("C55", idx)};

        return new Matrix(tmpArr, 5);
    }

    /**
     * Prints the data associated to one cross along with its first segment, cluster and hit.
     * @param cross the cross to be printed
     */
    @SuppressWarnings("unused")
    public static void printSampleCross(Cross cross) {
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
    @SuppressWarnings("unused")
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
