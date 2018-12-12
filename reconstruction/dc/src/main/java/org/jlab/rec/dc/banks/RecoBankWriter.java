package org.jlab.rec.dc.banks;

import java.util.ArrayList;
import java.util.List;

import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.io.hipo.HipoDataEvent;
import org.jlab.jnp.hipo.data.HipoEvent;

import org.jlab.rec.dc.hit.Hit;
import org.jlab.rec.dc.hit.FittedHit;
import org.jlab.rec.dc.cluster.FittedCluster;
import org.jlab.rec.dc.segment.Segment;
import org.jlab.rec.dc.cross.Cross;
import org.jlab.rec.dc.cross.CrossList;
import org.jlab.rec.dc.track.Track;

import trackfitter.fitter.utilities.*;

/**
 * A class to fill the reconstructed DC banks.
 *
 * @author ziegler
 */
public class RecoBankWriter {

    /**
     * Transforms a list hits into one of fitted hits.
     * @param hits list of hits to be transformed
     * @return     list of fitted hits
     */
    public List<FittedHit> createRawHitList(List<Hit> hits) {

        List<FittedHit> fhits = new ArrayList<>();

        for (Hit hit : hits) {
            FittedHit fhit = new FittedHit(hit.get_Sector(),
                                           hit.get_Superlayer(),
                                           hit.get_Layer(),
                                           hit.get_Wire(),
                                           hit.get_TDC(),
                                           hit.get_Id());
            fhit.set_Id(hit.get_Id());
            fhit.set_DocaErr(hit.get_DocaErr());
            fhits.add(fhit);
        }
        return fhits;
    }

    /**
     * Updates the list of fitted hits with data from the clusters containing them.
     * @param fhits    the list of fitted hits
     * @param clusters the list of clusters
     */
    public void updateHitsListWithClusterInfo(List<FittedHit> fhits, List<FittedCluster> clusters) {

        for (int i = 0; i < clusters.size(); i++) {
            clusters.get(i).set_Id(i + 1);
            for (int j = 0; j < clusters.get(i).size(); j++) {
                clusters.get(i).get(j).set_AssociatedClusterID(clusters.get(i).get_Id());
                for (int k = 0; k < fhits.size(); k++) {
                    if (fhits.get(k).get_Id() == clusters.get(i).get(j).get_Id()) {
                        fhits.remove(k);
                        fhits.add(clusters.get(i).get(j));
                    }
                }
            }
        }
    }

    /**
     * Writes a list of hits into the EvioEvent's HB bank.
     * @param event   the EvioEvent
     * @param hitlist the list of hits
     * @param TB      boolean set to 1 if time-based and 0 otherwise.
     * @return        hits bank
     */
    private DataBank fillHitsBank(DataEvent event, List<FittedHit> hitlist, boolean TB) {

        if (TB && event.hasBank("TimeBasedTrkg::TBHits")) { // For second pass tracking
            ((HipoDataEvent) event).getHipoEvent().removeGroup("TimeBasedTrkg::TBHits");
        }

        DataBank bank;
        if (!TB) bank = event.createBank("HitBasedTrkg::HBHits",  hitlist.size());
        else     bank = event.createBank("TimeBasedTrkg::TBHits", hitlist.size());

        for (int i = 0; i < hitlist.size(); i++) {
            if (hitlist.get(i).get_Id() == -1) continue;
            if (TB && hitlist.get(i).get_TrkResid() == 999) {
                hitlist.get(i).set_AssociatedTBTrackID(-1);
            }

            bank.setShort("id",         i, (short) hitlist.get(i).get_Id());
            bank.setByte ("superlayer", i, (byte)  hitlist.get(i).get_Superlayer());
            bank.setByte ("layer",      i, (byte)  hitlist.get(i).get_Layer());
            bank.setByte ("sector",     i, (byte)  hitlist.get(i).get_Sector());
            bank.setShort("wire",       i, (short) hitlist.get(i).get_Wire());
            bank.setFloat("X",          i, (float) hitlist.get(i).get_X());
            bank.setFloat("Z",          i, (float) hitlist.get(i).get_Z());
            bank.setByte ("LR",         i, (byte)  hitlist.get(i).get_LeftRightAmb());

            bank.setFloat("docaError",  i, (float) hitlist.get(i).get_DocaErr());
            bank.setFloat("trkDoca",    i, (float) hitlist.get(i).get_ClusFitDoca());
            bank.setShort("clusterID",  i, (short) hitlist.get(i).get_AssociatedClusterID());
            bank.setByte ("trkID",      i, (byte)  hitlist.get(i).get_AssociatedHBTrackID());

            bank.setInt  ("TDC",          i,         hitlist.get(i).get_TDC());
            bank.setFloat("B",            i, (float) hitlist.get(i).getB());
            bank.setFloat("TProp",        i, (float) hitlist.get(i).getTProp());
            bank.setFloat("TFlight",      i, (float) hitlist.get(i).getTFlight());

            // bank.setFloat("time",         i, (float) hitlist.get(i).get_Time());

            if (!TB) {
                bank.setShort("status", i, (short) 0);
                bank.setFloat("LocX",   i, (float) hitlist.get(i).get_lX());
                bank.setFloat("LocY",   i, (float) hitlist.get(i).get_lY());

                if (hitlist.get(i).get_AssociatedHBTrackID() > -1 && !event.hasBank("MC::Particle")) {
                    bank.setFloat("TProp", i, (float) hitlist.get(i).getSignalPropagTimeAlongWire());
                    bank.setFloat("TFlight", i, (float) hitlist.get(i).getSignalTimeOfFlight());
                }
                continue;
            }

            // else
            bank.setShort("status", i, (short) hitlist.get(i).get_QualityFac());

            // Checks the existing schema to fill the time
            if (bank.getDescriptor().hasEntry("time")) {
                bank.setFloat("time", i, (float) (hitlist.get(i).get_Time() -
                                                  hitlist.get(i).get_DeltaTimeBeta()));
            }
            if (bank.getDescriptor().hasEntry("beta")) {
                bank.setFloat("beta", i, (float) hitlist.get(i).get_Beta());
            }
            if (bank.getDescriptor().hasEntry("tBeta")) {
                bank.setFloat("tBeta", i, (float) hitlist.get(i).get_DeltaTimeBeta());
            }
            if (bank.getDescriptor().hasEntry("fitResidual")) {
                bank.setFloat("fitResidual", i, (float) hitlist.get(i).get_TrkResid());
            }

            bank.setFloat("doca",         i, (float) hitlist.get(i).get_Doca());
            bank.setFloat("timeResidual", i, (float) hitlist.get(i).get_TimeResidual());

            bank.setFloat("T0",           i, (float) hitlist.get(i).getT0());
            bank.setFloat("TStart",       i, (float) hitlist.get(i).getTStart());

            if (hitlist.get(i).get_AssociatedTBTrackID() > -1 && !event.hasBank("MC::Particle")) {
                if (hitlist.get(i).getSignalPropagTimeAlongWire() != 0 &&
                        hitlist.get(i).get_AssociatedTBTrackID() >= 1) {
                    bank.setFloat("TProp", i, (float) hitlist.get(i).getSignalPropagTimeAlongWire());
                }

                if (hitlist.get(i).getSignalTimeOfFlight() != 0 &&
                        hitlist.get(i).get_AssociatedTBTrackID() >= 1) {
                    bank.setFloat("TFlight", i, (float) hitlist.get(i).getSignalTimeOfFlight());
                }
            }
        }

        return bank;
    }

    /**
     * Writes a list of clusters into the EvioEvent's bank.
     * @param event    the EvioEvent
     * @param cluslist the list of clusters
     * @param TB       boolean set to 1 if time-based and 0 otherwise.
     * @return         clusters bank
     */
    private DataBank fillClustersBank(DataEvent event, List<FittedCluster> cluslist, boolean TB) {

        if (TB && event.hasBank("TimeBasedTrkg::TBClusters")) {
            // For second pass tracking
            ((HipoDataEvent) event).getHipoEvent().removeGroup("TimeBasedTrkg::TBClusters");
        }

        DataBank bank;
        if (!TB) bank = event.createBank("HitBasedTrkg::HBClusters",  cluslist.size());
        else     bank = event.createBank("TimeBasedTrkg::TBClusters", cluslist.size());

        for (int i = 0; i < cluslist.size(); i++) {
            if (cluslist.get(i).get_Id() == -1) continue;
            List<Integer> hitIdxArray = new ArrayList<Integer>();

            double chi2 = 0;

            bank.setShort("id",         i, (short) cluslist.get(i).get_Id());
            bank.setByte ("superlayer", i, (byte)  cluslist.get(i).get_Superlayer());
            bank.setByte ("sector",     i, (byte)  cluslist.get(i).get_Sector());
            bank.setFloat("avgWire",    i, (float) cluslist.get(i).getAvgwire());
            bank.setByte ("size",       i, (byte)  cluslist.get(i).size());

            if (!TB) {
                if (cluslist.get(i).size() < 6) bank.setShort("status", i, (short) 1);
                else                            bank.setShort("status", i, (short) 0);
            }

            double fitSlope  = cluslist.get(i).get_clusterLineFitSlope();
            double fitInterc = cluslist.get(i).get_clusterLineFitIntercept();

            bank.setFloat("fitSlope",     i, (float) fitSlope);
            bank.setFloat("fitSlopeErr",  i, (float) cluslist.get(i).get_clusterLineFitSlopeErr());
            bank.setFloat("fitInterc",    i, (float) fitInterc);
            bank.setFloat("fitIntercErr", i, (float) cluslist.get(i).get_clusterLineFitInterceptErr());

            for (int j = 0; j < cluslist.get(i).size() && j < 12; j++) {
                hitIdxArray.add(cluslist.get(i).get(j).get_Id());

                // TODO: See further into these lines
                // Math.sqrt(12.) = 3.4641016151377544
                // double residual = cluslist.get(i).get(j).get_ClusFitDoca() /
                //                   (cluslist.get(i).get(j).get_CellSize() / 3.4641016151377544);
                // chi2 += residual * residual;
            }
            // bank.setFloat("fitChisqProb", i,
            //               (float) ProbChi2perNDF.prob(chi2, cluslist.get(i).size() - 2));
            bank.setFloat("fitChisqProb", i, (float) cluslist.get(i).get_fitProb());

            for (int j = 1; j < hitIdxArray.size() + 1; j++) {
                bank.setShort("Hit" + j + "_ID", i, (short) hitIdxArray.get(j - 1).intValue());
            }
        }

        return bank;
    }

    /**
     * Writes a list of segments into the EvioEvent's bank.
     * @param event   the EvioEvent
     * @param seglist the list of segments
     * @param TB      boolean set to 1 if time-based and 0 otherwise.
     * @return        segments bank
     */
    private DataBank fillSegmentsBank(DataEvent event, List<Segment> seglist, boolean TB) {

        if (TB && event.hasBank("TimeBasedTrkg::TBSegments")) { // For second pass tracking
            ((HipoDataEvent) event).getHipoEvent().removeGroup("TimeBasedTrkg::TBSegments");
        }

        DataBank bank;
        if (!TB) bank = event.createBank("HitBasedTrkg::HBSegments",  seglist.size());
        else     bank = event.createBank("TimeBasedTrkg::TBSegments", seglist.size());

        int[] hitIdxArray = new int[12];

        for (int i = 0; i < seglist.size(); i++) {
            if (seglist.get(i).get_Id() == -1) continue;

            for (int j = 0; j < hitIdxArray.length; j++) {
                hitIdxArray[j] = -1;
            }

            double chi2 = 0;

            bank.setShort("id",        i, (short) seglist.get(i).get_Id());
            bank.setByte("superlayer", i, (byte) seglist.get(i).get_Superlayer());
            bank.setByte("sector",     i, (byte) seglist.get(i).get_Sector());
            if (TB) bank.setShort("status", i, (short) seglist.get(i).get_Status());

            FittedCluster cls = seglist.get(i).get_fittedCluster();
            bank.setShort("Cluster_ID",    i, (short) cls.get_Id());

            bank.setFloat("avgWire",       i, (float) cls.getAvgwire());
            bank.setByte ("size",          i, (byte)  seglist.get(i).size());
            bank.setFloat("fitSlope",      i, (float) cls.get_clusterLineFitSlope());
            bank.setFloat("fitSlopeErr",   i, (float) cls.get_clusterLineFitSlopeErr());
            bank.setFloat("fitInterc",     i, (float) cls.get_clusterLineFitIntercept());
            bank.setFloat("fitIntercErr",  i, (float) cls.get_clusterLineFitInterceptErr());
            bank.setFloat("SegEndPoint1X", i, (float) seglist.get(i).get_SegmentEndPoints()[0]);
            bank.setFloat("SegEndPoint1Z", i, (float) seglist.get(i).get_SegmentEndPoints()[1]);
            bank.setFloat("SegEndPoint2X", i, (float) seglist.get(i).get_SegmentEndPoints()[2]);
            bank.setFloat("SegEndPoint2Z", i, (float) seglist.get(i).get_SegmentEndPoints()[3]);

            if (TB) {
                bank.setFloat("resiSum", i, (float) seglist.get(i).get_ResiSum());
                bank.setFloat("timeSum", i, (float) seglist.get(i).get_TimeSum());
            }

            for (int j = 0; j < seglist.get(i).size(); j++) {
                if (seglist.get(i).get_Id() == -1) continue;
                if (j < hitIdxArray.length) {
                    hitIdxArray[j] = seglist.get(i).get(j).get_Id();
                }
                // Math.sqrt(12.) = 3.4641016151377544
                double residual = seglist.get(i).get(j).get_ClusFitDoca() /
                                  (seglist.get(i).get(j).get_CellSize() / 3.4641016151377544);
                chi2 += residual * residual;
            }
            bank.setFloat("fitChisqProb", i,
                          (float) ProbChi2perNDF.prob(chi2, seglist.get(i).size() - 2));

            for (int j = 0; j < hitIdxArray.length; j++) {
                String hitStrg = "Hit";
                hitStrg += (j + 1);
                hitStrg += "_ID";
                bank.setShort(hitStrg, i, (short) hitIdxArray[j]);
            }
        }

        return bank;
    }

    /**
     * Writes a list of crosses into the EvioEvent's bank.
     * @param event     the EvioEvent
     * @param crossList the list of crosses
     * @param TB        boolean set to 1 if time-based and 0 otherwise.
     * @return          crosses bank
     */
    private DataBank fillCrossesBank(DataEvent event, List<Cross> crosslist, boolean TB) {

        if (TB && event.hasBank("TimeBasedTrkg::TBCrosses")) { // For second pass tracking
            ((HipoDataEvent) event).getHipoEvent().removeGroup("TimeBasedTrkg::TBCrosses");
        }
        int banksize = 0;
        for (Cross cross : crosslist) {
            if (cross.get_Id() != -1) banksize++;
        }

        DataBank bank;
        if (!TB) bank = event.createBank("HitBasedTrkg::HBCrosses",  banksize);
        else     bank = event.createBank("TimeBasedTrkg::TBCrosses", banksize);

        int idx = 0;
        for (Cross cross : crosslist) {
            if (cross.get_Id() == -1) continue;
            bank.setShort("id",          idx, (short) cross.get_Id());
            bank.setByte ("sector",      idx, (byte)  cross.get_Sector());
            bank.setByte ("region",      idx, (byte)  cross.get_Region());
            bank.setFloat("x",           idx, (float) cross.get_Point().x());
            bank.setFloat("y",           idx, (float) cross.get_Point().y());
            bank.setFloat("z",           idx, (float) cross.get_Point().z());
            bank.setFloat("err_x",       idx, (float) cross.get_PointErr().x());
            bank.setFloat("err_y",       idx, (float) cross.get_PointErr().y());
            bank.setFloat("err_z",       idx, (float) cross.get_PointErr().z());
            bank.setFloat("ux",          idx, (float) cross.get_Dir().x());
            bank.setFloat("uy",          idx, (float) cross.get_Dir().y());
            bank.setFloat("uz",          idx, (float) cross.get_Dir().z());
            bank.setFloat("err_ux",      idx, (float) cross.get_DirErr().x());
            bank.setFloat("err_uy",      idx, (float) cross.get_DirErr().y());
            bank.setFloat("err_uz",      idx, (float) cross.get_DirErr().z());
            bank.setShort("Segment1_ID", idx, (short) cross.get_Segment1().get_Id());
            bank.setShort("Segment2_ID", idx, (short) cross.get_Segment2().get_Id());

            if (!TB) bank.setShort("status", idx, (short) 0);
            else     bank.setShort("status", idx, (short) (cross.get_Segment1().get_Status() +
                                                           cross.get_Segment2().get_Status()));
            idx++;
        }

        return bank;
    }

    /**
     * Writes a list of tracks into the EvioEvent's bank.
     * @param event    the EvioEvent
     * @param candlist the list of tracks
     * @param TB       boolean set to 1 if time-based and 0 otherwise.
     * @return         tracks bank
     */
    private DataBank fillTracksBank(DataEvent event, List<Track> candlist, boolean TB) {

        if (TB && event.hasBank("TimeBasedTrkg::TBTracks")) { // For second pass tracking
            ((HipoDataEvent) event).getHipoEvent().removeGroup("TimeBasedTrkg::TBTracks");
        }

        DataBank bank;
        if (!TB) bank = event.createBank("HitBasedTrkg::HBTracks",  candlist.size());
        else     bank = event.createBank("TimeBasedTrkg::TBTracks", candlist.size());

        for (int i = 0; i < candlist.size(); i++) {
            bank.setShort("id",     i, (short) candlist.get(i).get_Id());
            bank.setByte ("sector", i, (byte)  candlist.get(i).get_Sector());
            bank.setByte ("q",      i, (byte)  candlist.get(i).get_Q());
            bank.setShort("status", i, (short) (100 + candlist.get(i).get_Status() * 10 +
                                                candlist.get(i).get_MissingSuperlayer()));

            if (candlist.get(i).get_PreRegion1CrossPoint() != null) {
                bank.setFloat("c1_x",  i, (float) candlist.get(i).get_PreRegion1CrossPoint().x());
                bank.setFloat("c1_y",  i, (float) candlist.get(i).get_PreRegion1CrossPoint().y());
                bank.setFloat("c1_z",  i, (float) candlist.get(i).get_PreRegion1CrossPoint().z());
                bank.setFloat("c1_ux", i, (float) candlist.get(i).get_PreRegion1CrossDir().x());
                bank.setFloat("c1_uy", i, (float) candlist.get(i).get_PreRegion1CrossDir().y());
                bank.setFloat("c1_uz", i, (float) candlist.get(i).get_PreRegion1CrossDir().z());
            }
            if (candlist.get(i).get_PostRegion3CrossPoint() != null) {
                bank.setFloat("c3_x",  i, (float) candlist.get(i).get_PostRegion3CrossPoint().x());
                bank.setFloat("c3_y",  i, (float) candlist.get(i).get_PostRegion3CrossPoint().y());
                bank.setFloat("c3_z",  i, (float) candlist.get(i).get_PostRegion3CrossPoint().z());
                bank.setFloat("c3_ux", i, (float) candlist.get(i).get_PostRegion3CrossDir().x());
                bank.setFloat("c3_uy", i, (float) candlist.get(i).get_PostRegion3CrossDir().y());
                bank.setFloat("c3_uz", i, (float) candlist.get(i).get_PostRegion3CrossDir().z());
            }
            if (candlist.get(i).get_Region1TrackX() != null) {
                bank.setFloat("t1_x",  i, (float) candlist.get(i).get_Region1TrackX().x());
                bank.setFloat("t1_y",  i, (float) candlist.get(i).get_Region1TrackX().y());
                bank.setFloat("t1_z",  i, (float) candlist.get(i).get_Region1TrackX().z());
                bank.setFloat("t1_px", i, (float) candlist.get(i).get_Region1TrackP().x());
                bank.setFloat("t1_py", i, (float) candlist.get(i).get_Region1TrackP().y());
                bank.setFloat("t1_pz", i, (float) candlist.get(i).get_Region1TrackP().z());
            }

            bank.setFloat("pathlength", i, (float) candlist.get(i).get_TotPathLen());
            bank.setFloat("Vtx0_x",     i, (float) candlist.get(i).get_Vtx0().x());
            bank.setFloat("Vtx0_y",     i, (float) candlist.get(i).get_Vtx0().y());
            bank.setFloat("Vtx0_z",     i, (float) candlist.get(i).get_Vtx0().z());
            bank.setFloat("p0_x",       i, (float) candlist.get(i).get_pAtOrig().x());
            bank.setFloat("p0_y",       i, (float) candlist.get(i).get_pAtOrig().y());
            bank.setFloat("p0_z",       i, (float) candlist.get(i).get_pAtOrig().z());
            bank.setFloat("chi2",       i, (float) candlist.get(i).get_FitChi2());
            bank.setShort("ndf",        i, (short) candlist.get(i).get_FitNDF());

            if (!TB) {
                bank.setShort("Cross1_ID",  i, (short) candlist.get(i).get(0).get_Id());
                bank.setShort("Cross2_ID",  i, (short) candlist.get(i).get(1).get_Id());
                bank.setShort("Cross3_ID",  i, (short) candlist.get(i).get(2).get_Id());
                continue;
            }

            // else
            if (candlist.get(i).size() == 3) {
                bank.setShort("Cross1_ID", i, (short) candlist.get(i).get(0).get_Id());
                bank.setShort("Cross2_ID", i, (short) candlist.get(i).get(1).get_Id());
                bank.setShort("Cross3_ID", i, (short) candlist.get(i).get(2).get_Id());
            }
            else if (candlist.get(i).size() == 2) {
                bank.setShort("Cross1_ID", i, (short) candlist.get(i).get(0).get_Id());
                bank.setShort("Cross2_ID", i, (short) candlist.get(i).get(1).get_Id());
                bank.setShort("Cross3_ID", i, (short) -1);
            }
            else if (candlist.get(i).size() == 1) {
                bank.setShort("Cross1_ID", i, (short) candlist.get(i).get(0).get_Id());
                bank.setShort("Cross2_ID", i, (short) -1);
                bank.setShort("Cross3_ID", i, (short) -1);
            }
        }
        return bank;
    }

    /**
     * writes the covariance matrix from HB fits to be used for starting the Time-based tracking.
     * @param event    hipo event
     * @param candlist list of tracks
     * @return         covariance matrix
     */
    private DataBank fillTrackCovMatBank(DataEvent event, List<Track> candlist) {

        if (event == null) return null;
        DataBank bank = event.createBank("TimeBasedTrkg::TBCovMat", candlist.size());

        for (int i = 0; i < candlist.size(); i++) {
            bank.setShort("id", i, (short) candlist.get(i).get_Id());
            if (candlist.get(i).get_CovMat() == null) continue;
            bank.setFloat("C11", i, (float) candlist.get(i).get_CovMat().get(0, 0));
            bank.setFloat("C12", i, (float) candlist.get(i).get_CovMat().get(0, 1));
            bank.setFloat("C13", i, (float) candlist.get(i).get_CovMat().get(0, 2));
            bank.setFloat("C14", i, (float) candlist.get(i).get_CovMat().get(0, 3));
            bank.setFloat("C15", i, (float) candlist.get(i).get_CovMat().get(0, 4));
            bank.setFloat("C21", i, (float) candlist.get(i).get_CovMat().get(1, 0));
            bank.setFloat("C22", i, (float) candlist.get(i).get_CovMat().get(1, 1));
            bank.setFloat("C23", i, (float) candlist.get(i).get_CovMat().get(1, 2));
            bank.setFloat("C24", i, (float) candlist.get(i).get_CovMat().get(1, 3));
            bank.setFloat("C25", i, (float) candlist.get(i).get_CovMat().get(1, 4));
            bank.setFloat("C31", i, (float) candlist.get(i).get_CovMat().get(2, 0));
            bank.setFloat("C32", i, (float) candlist.get(i).get_CovMat().get(2, 1));
            bank.setFloat("C33", i, (float) candlist.get(i).get_CovMat().get(2, 2));
            bank.setFloat("C34", i, (float) candlist.get(i).get_CovMat().get(2, 3));
            bank.setFloat("C35", i, (float) candlist.get(i).get_CovMat().get(2, 4));
            bank.setFloat("C41", i, (float) candlist.get(i).get_CovMat().get(3, 0));
            bank.setFloat("C42", i, (float) candlist.get(i).get_CovMat().get(3, 1));
            bank.setFloat("C43", i, (float) candlist.get(i).get_CovMat().get(3, 2));
            bank.setFloat("C44", i, (float) candlist.get(i).get_CovMat().get(3, 3));
            bank.setFloat("C45", i, (float) candlist.get(i).get_CovMat().get(3, 4));
            bank.setFloat("C51", i, (float) candlist.get(i).get_CovMat().get(4, 0));
            bank.setFloat("C52", i, (float) candlist.get(i).get_CovMat().get(4, 1));
            bank.setFloat("C53", i, (float) candlist.get(i).get_CovMat().get(4, 2));
            bank.setFloat("C54", i, (float) candlist.get(i).get_CovMat().get(4, 3));
            bank.setFloat("C55", i, (float) candlist.get(i).get_CovMat().get(4, 4));
        }

        return bank;
    }

    /**
     * Writes a list of trajectories with data pulled from a list of tracks.
     * @param event  the EvioEvent
     * @param tracks the list of tracks
     * @return       trajectories bank
     */
    private DataBank fillTrajectoryBank(DataEvent event, List<Track> tracks) {
        DataBank bank = event.createBank("TimeBasedTrkg::Trajectory", tracks.size() * 21);

        int idx = 0;
        for (Track track : tracks) {
            if (track == null)            continue;
            if (track.trajectory == null) continue;

            for (int j = 0; j < track.trajectory.size(); j++) {
                if (track.trajectory.get(j).getDetName().startsWith("DC") && (j - 6) % 6 != 0) {
                    continue; // save the last layer in a superlayer
                }

                bank.setShort("did", idx, (short) track.trajectory.get(j).getDetId());
                bank.setShort("tid", idx, (short) track.get_Id());
                bank.setFloat("x",   idx, (float) track.trajectory.get(j).getX());
                bank.setFloat("y",   idx, (float) track.trajectory.get(j).getY());
                bank.setFloat("z",   idx, (float) track.trajectory.get(j).getZ());
                bank.setFloat("tx",  idx, (float) ((float) track.trajectory.get(j).getpX() /
                                                   track.get_P()));
                bank.setFloat("ty",  idx, (float) ((float) track.trajectory.get(j).getpY() /
                                                   track.get_P()));
                bank.setFloat("tz", idx,  (float) ((float) track.trajectory.get(j).getpZ() /
                                                   track.get_P()));
                bank.setFloat("B", idx,   (float) track.trajectory.get(j).getiBdl());
                bank.setFloat("L", idx,   (float) track.trajectory.get(j).getPathLen());

                idx++;
            }
        }

        return bank;
    }

    /**
     * Writes the hits, clusters, segments, crosses and track candidates into the EvioEvent's
     * Hit-based or Time-based bank.
     * @param event    hipo event
     * @param rbw      RecoBankWriter's instance where everything is written to
     * @param fhits    list of hits
     * @param clusters list of clusters
     * @param segments list of segments
     * @param crosses  list of crosses
     * @param trkcands list of tracks
     * @param TB       boolean set to 1 if time-based and 0 otherwise.
     */
    private void fillAllBanks(DataEvent event, RecoBankWriter rbw,
                             List<FittedHit>     fhits,
                             List<FittedCluster> clusters,
                             List<Segment>       segments,
                             List<Cross>         crosses,
                             List<Track>         trkcands,
                             boolean             TB) {

        if (event == null) return;
        if (fhits != null) event.appendBanks(rbw.fillHitsBank(event, fhits, TB));
        else return;

        if (clusters != null) event.appendBanks(rbw.fillClustersBank(event, clusters, TB));
        else return;

        if (segments != null) event.appendBanks(rbw.fillSegmentsBank(event, segments, TB));
        else return;

        if (crosses != null) event.appendBanks(rbw.fillCrossesBank(event, crosses, TB));
        else return;

        if (trkcands != null) {
            event.appendBanks(rbw.fillTracksBank(event, trkcands, TB));
            if (!TB) event.appendBanks(rbw.fillTrackCovMatBank(event, trkcands));
            else     event.appendBanks(rbw.fillTrajectoryBank(event, trkcands));
        }
        return;
    }

    public void fillAllHBBanks(DataEvent event, RecoBankWriter rbw,
                               List<FittedHit>     fhits,
                               List<FittedCluster> clusters,
                               List<Segment>       segments,
                               List<Cross>         crosses,
                               List<Track>         trkcands) {
        fillAllBanks(event, rbw, fhits, clusters, segments, crosses, trkcands, false);
        return;
    }

    public void fillHBTracksBanks(DataEvent event, RecoBankWriter rbw, List<Track> trkcands) {
        if (event.hasBank("HitBasedTrkg::HBTracks"))
            ((HipoDataEvent) event).getHipoEvent().removeGroup("HitBasedTrkg::HBTracks");
        if (event.hasBank("TimeBasedTrkg::TBCovMat"))
            ((HipoDataEvent) event).getHipoEvent().removeGroup("TimeBasedTrkg::TBCovMat");

        if (trkcands != null) {
            event.appendBanks(rbw.fillTracksBank(event, trkcands, false));
            event.appendBanks(rbw.fillTrackCovMatBank(event, trkcands));
        }

        return;
    }

    public void fillAllTBBanks(DataEvent event, RecoBankWriter rbw,
                               List<FittedHit>     fhits,
                               List<FittedCluster> clusters,
                               List<Segment>       segments,
                               List<Cross>         crosses,
                               List<Track>         trkcands) {
        fillAllBanks(event, rbw, fhits, clusters, segments, crosses, trkcands, true);
        return;
    }
}
