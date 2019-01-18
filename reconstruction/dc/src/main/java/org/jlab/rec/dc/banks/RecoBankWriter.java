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
import org.jlab.rec.dc.trajectory.StateVec;
import org.jlab.rec.dc.track.Track;

import trackfitter.fitter.utilities.*;

/**
 * A class to fill the reconstructed DC banks
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
            FittedHit fhit = new FittedHit(hit.get_Sector(), hit.get_Superlayer(),
                                           hit.get_Layer(),  hit.get_Wire(),
                                           hit.get_TDC(),    hit.get_Id());
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
     * @param hitList the list of hits
     * @param TB      boolean set to 1 if time-based and 0 otherwise.
     * @return        hits bank
     */
    private DataBank fillHitsBank(DataEvent event, List<FittedHit> hitList, boolean TB) {
        if (TB && event.hasBank("TimeBasedTrkg::TBHits")) { // For second pass tracking
            ((HipoDataEvent) event).getHipoEvent().removeGroup("TimeBasedTrkg::TBHits");
        }
        if (!TB && event.hasBank("HitBasedTrkg::HBHits")) {
            ((HipoDataEvent) event).getHipoEvent().removeGroup("HitBasedTrkg::HBHits");
        }

        DataBank bank;
        if (!TB) bank = event.createBank("HitBasedTrkg::HBHits",  hitList.size());
        else     bank = event.createBank("TimeBasedTrkg::TBHits", hitList.size());

        for (int i = 0; i < hitList.size(); i++) {
            if (hitList.get(i).get_Id() == -1) continue;
            if (TB && hitList.get(i).get_TrkResid() == 999) {
                hitList.get(i).set_AssociatedTBTrackID(-1);
            }

            bank.setShort("id",         i, (short) hitList.get(i).get_Id());
            bank.setByte ("superlayer", i, (byte)  hitList.get(i).get_Superlayer());
            bank.setByte ("layer",      i, (byte)  hitList.get(i).get_Layer());
            bank.setByte ("sector",     i, (byte)  hitList.get(i).get_Sector());
            bank.setShort("wire",       i, (short) hitList.get(i).get_Wire());
            bank.setFloat("X",          i, (float) hitList.get(i).get_X());
            bank.setFloat("Z",          i, (float) hitList.get(i).get_Z());
            bank.setByte ("LR",         i, (byte)  hitList.get(i).get_LeftRightAmb());

            bank.setFloat("docaError",  i, (float) hitList.get(i).get_DocaErr());
            bank.setFloat("trkDoca",    i, (float) hitList.get(i).get_ClusFitDoca());
            bank.setShort("clusterID",  i, (short) hitList.get(i).get_AssociatedClusterID());
            bank.setByte ("trkID",      i, (byte)  hitList.get(i).get_AssociatedHBTrackID());

            bank.setInt  ("TDC",          i,         hitList.get(i).get_TDC());
            bank.setFloat("B",            i, (float) hitList.get(i).getB());
            bank.setFloat("TProp",        i, (float) hitList.get(i).getTProp());
            bank.setFloat("TFlight",      i, (float) hitList.get(i).getTFlight());

            if (!TB) {
                bank.setFloat("cellSize", i, (float) hitList.get(i).get_CellSize());
                bank.setFloat("XMP", i, (float) hitList.get(i).get_XMP());
                bank.setFloat("residual", i, (float) hitList.get(i).get_Residual());
                bank.setFloat("timeResidual", i, (float) hitList.get(i).get_TimeResidual());
                bank.setFloat("qualityFac", i, (float) hitList.get(i).get_QualityFac());
                bank.setByte ("trkStatus", i, (byte) hitList.get(i).get_TrkgStatus());
                bank.setFloat("clusFitDoca", i, (float) hitList.get(i).get_ClusFitDoca());
                bank.setFloat("trkFitDoca", i, (float) hitList.get(i).get_TrkFitDoca());
                bank.setFloat("timeToDistance", i, (float) hitList.get(i).get_TimeToDistance());
                bank.setFloat("beta", i, (float) hitList.get(i).get_Beta());
                bank.setFloat("doca", i, (float) hitList.get(i).get_Doca());
                bank.setByte ("lr", i, (byte) hitList.get(i)._lr);
                // bank.setFloat("crossDirIntersWireX", i, (float) hitList.get(i).getCrossDirIntersWire().x());
                // bank.setFloat("crossDirIntersWireY", i, (float) hitList.get(i).getCrossDirIntersWire().y());
                // bank.setFloat("crossDirIntersWireZ", i, (float) hitList.get(i).getCrossDirIntersWire().z());
                bank.setFloat("signalPropagAlongWire", i, (float) hitList.get(i).getSignalPropagAlongWire());
                bank.setFloat("signalTimeOfFlight", i, (float) hitList.get(i).getSignalTimeOfFlight());
                bank.setFloat("t0", i, (float) hitList.get(i).getT0());
                bank.setFloat("tStart", i, (float) hitList.get(i).getTStart());
                bank.setFloat("time", i, (float) hitList.get(i).get_Time());
                bank.setByte ("outOfTimeFlag", i, (byte) (hitList.get(i).get_OutOfTimeFlag() ? 1 : 0));
                bank.setFloat("wireLength", i, (float) hitList.get(i).get_WireLength());
                bank.setFloat("wireMaxSag", i, (float) hitList.get(i).get_WireMaxSag());
                bank.setFloat("trkResid", i, (float) hitList.get(i).get_TrkResid());
                bank.setFloat("deltaTimeBeta", i, (float) hitList.get(i).get_DeltaTimeBeta());
            }

            if (!TB) {
                bank.setShort("status", i, (short) 0);
                bank.setFloat("LocX",   i, (float) hitList.get(i).get_lX());
                bank.setFloat("LocY",   i, (float) hitList.get(i).get_lY());

                if (hitList.get(i).get_AssociatedHBTrackID() > -1 && !event.hasBank("MC::Particle")) {
                    bank.setFloat("TProp", i, (float) hitList.get(i).getSignalPropagTimeAlongWire());
                    bank.setFloat("TFlight", i, (float) hitList.get(i).getSignalTimeOfFlight());
                }
                continue;
            }

            // else
            bank.setShort("status", i, (short) hitList.get(i).get_QualityFac());

            // Checks the existing schema to fill the time
            if (bank.getDescriptor().hasEntry("time")) {
                bank.setFloat("time", i, (float) (hitList.get(i).get_Time() -
                                                  hitList.get(i).get_DeltaTimeBeta()));
            }
            if (bank.getDescriptor().hasEntry("beta")) {
                bank.setFloat("beta", i, (float) hitList.get(i).get_Beta());
            }
            if (bank.getDescriptor().hasEntry("tBeta")) {
                bank.setFloat("tBeta", i, (float) hitList.get(i).get_DeltaTimeBeta());
            }
            if (bank.getDescriptor().hasEntry("fitResidual")) {
                bank.setFloat("fitResidual", i, (float) hitList.get(i).get_TrkResid());
            }

            bank.setFloat("doca",         i, (float) hitList.get(i).get_Doca());
            bank.setFloat("timeResidual", i, (float) hitList.get(i).get_TimeResidual());

            bank.setFloat("T0",           i, (float) hitList.get(i).getT0());
            bank.setFloat("TStart",       i, (float) hitList.get(i).getTStart());

            if (hitList.get(i).get_AssociatedTBTrackID() > -1 && !event.hasBank("MC::Particle")) {
                if (hitList.get(i).getSignalPropagTimeAlongWire() != 0 &&
                        hitList.get(i).get_AssociatedTBTrackID() >= 1) {
                    bank.setFloat("TProp", i, (float) hitList.get(i).getSignalPropagTimeAlongWire());
                }

                if (hitList.get(i).getSignalTimeOfFlight() != 0 &&
                        hitList.get(i).get_AssociatedTBTrackID() >= 1) {
                    bank.setFloat("TFlight", i, (float) hitList.get(i).getSignalTimeOfFlight());
                }
            }
        }

        return bank;
    }

    /**
     * Writes a list of clusters into the EvioEvent's bank.
     * @param event    the EvioEvent
     * @param clusList the list of clusters
     * @param TB       boolean set to 1 if time-based and 0 otherwise.
     * @return         clusters bank
     */
    private DataBank fillClustersBank(DataEvent event, List<FittedCluster> clusList, boolean TB,
                                      boolean DCHB2) {
        if (TB && event.hasBank("TimeBasedTrkg::TBClusters")) {
            // For second pass tracking
            ((HipoDataEvent) event).getHipoEvent().removeGroup("TimeBasedTrkg::TBClusters");
        }
        if (!TB && event.hasBank("HitBasedTrkg::HBClusters")) {
            ((HipoDataEvent) event).getHipoEvent().removeGroup("HitBasedTrkg::HBClusters");
        }

        DataBank bank;
        if (!TB) bank = event.createBank("HitBasedTrkg::HBClusters",  clusList.size());
        else     bank = event.createBank("TimeBasedTrkg::TBClusters", clusList.size());

        for (int i = 0; i < clusList.size(); i++) {
            if (clusList.get(i).get_Id() == -1) continue;
            List<Integer> hitIdxArray = new ArrayList<Integer>();

            double chi2 = 0;

            bank.setShort("id",         i, (short) clusList.get(i).get_Id());
            bank.setByte ("superlayer", i, (byte)  clusList.get(i).get_Superlayer());
            bank.setByte ("sector",     i, (byte)  clusList.get(i).get_Sector());
            bank.setFloat("avgWire",    i, (float) clusList.get(i).getAvgwire());
            bank.setByte ("size",       i, (byte)  clusList.get(i).size());

            if (!TB) {
                if (clusList.get(i).size() < 6) bank.setShort("status", i, (short) 1);
                else                            bank.setShort("status", i, (short) 0);
            }

            bank.setFloat("fitSlope",     i, (float) clusList.get(i).get_clusterLineFitSlope());
            bank.setFloat("fitSlopeErr",  i, (float) clusList.get(i).get_clusterLineFitSlopeErr());
            bank.setFloat("fitInterc",    i, (float) clusList.get(i).get_clusterLineFitIntercept());
            bank.setFloat("fitIntercErr", i, (float) clusList.get(i).get_clusterLineFitInterceptErr());

            for (int j = 0; j < clusList.get(i).size() && j < 12; j++) {
                hitIdxArray.add(clusList.get(i).get(j).get_Id());

                if (DCHB2) {
                    // Math.sqrt(12.) = 3.4641016151377544
                    double residual = clusList.get(i).get(j).get_ClusFitDoca() /
                                      (clusList.get(i).get(j).get_CellSize() / 3.4641016151377544);
                    chi2 += residual * residual;
                }
            }
            if (DCHB2) {
                bank.setFloat("fitChisqProb", i,
                              (float) ProbChi2perNDF.prob(chi2, clusList.get(i).size() - 2));
            }
            else {
                bank.setFloat("fitChisqProb", i, (float) clusList.get(i).get_fitProb());
            }
            for (int j = 1; j < hitIdxArray.size() + 1; j++) {
                bank.setShort("Hit" + j + "_ID", i, (short) hitIdxArray.get(j - 1).intValue());
            }

            if (!TB) {
                bank.setFloat("fitSlIntCov",  i, (float) clusList.get(i).get_clusterLineFitSlIntCov());
                bank.setFloat("chisqProb",    i, (float) clusList.get(i).get_Chisq());
                bank.setFloat("clusLine1X", i, (float) clusList.get(i).get_clusLine().origin().x());
                bank.setFloat("clusLine1Y", i, (float) clusList.get(i).get_clusLine().origin().y());
                bank.setFloat("clusLine1Z", i, (float) clusList.get(i).get_clusLine().origin().z());
                bank.setFloat("clusLine2X", i, (float) clusList.get(i).get_clusLine().end().x());
                bank.setFloat("clusLine2Y", i, (float) clusList.get(i).get_clusLine().end().y());
                bank.setFloat("clusLine2Z", i, (float) clusList.get(i).get_clusLine().end().z());
                bank.setFloat("clusLineErr1X", i, (float) clusList.get(i).get_clusLineErr().origin().x());
                bank.setFloat("clusLineErr1Y", i, (float) clusList.get(i).get_clusLineErr().origin().y());
                bank.setFloat("clusLineErr1Z", i, (float) clusList.get(i).get_clusLineErr().origin().z());
                bank.setFloat("clusLineErr2X", i, (float) clusList.get(i).get_clusLineErr().end().x());
                bank.setFloat("clusLineErr2Y", i, (float) clusList.get(i).get_clusLineErr().end().y());
                bank.setFloat("clusLineErr2Z", i, (float) clusList.get(i).get_clusLineErr().end().z());
                bank.setFloat("clusterLineFitSlopeMP", i,
                        (float) clusList.get(i).get_clusterLineFitSlopeMP());
                bank.setFloat("clusterLineFitSlopeErrMP", i,
                        (float) clusList.get(i).get_clusterLineFitSlopeErrMP());
                bank.setFloat("clusterLineFitInterceptMP", i,
                        (float) clusList.get(i).get_clusterLineFitInterceptMP());
                bank.setFloat("clusterLineFitInterceptErrMP", i,
                        (float) clusList.get(i).get_clusterLineFitInterceptErrMP());

                if (clusList.get(i).get_Status() != null) {
                    for (int c1 = 0; c1 < clusList.get(i).get_Status().length; ++c1) {
                        for (int c2 = 0; c2 < clusList.get(i).get_Status()[c1].length; ++c2) {
                            bank.setByte(("clusterStatus" + c1) + c2, i,
                                    (byte) clusList.get(i).get_Status()[c1][c2]);
                        }
                    }
                }
                else bank.setByte("clusterStatus00", i, (byte) -1);
            }
        }

        return bank;
    }

    /**
     * Writes a list of segments into the EvioEvent's bank.
     * @param event   the EvioEvent
     * @param segList the list of segments
     * @param TB      boolean set to 1 if time-based and 0 otherwise.
     * @return        segments bank
     */
    private DataBank fillSegmentsBank(DataEvent event, List<Segment> segList, boolean TB) {
        if (TB && event.hasBank("TimeBasedTrkg::TBSegments")) { // For second pass tracking
            ((HipoDataEvent) event).getHipoEvent().removeGroup("TimeBasedTrkg::TBSegments");
        }
        if (!TB && event.hasBank("HitBasedTrkg::HBSegments")) {
            ((HipoDataEvent) event).getHipoEvent().removeGroup("HitBasedTrkg::HBSegments");
        }

        DataBank bank;
        if (!TB) bank = event.createBank("HitBasedTrkg::HBSegments",  segList.size());
        else     bank = event.createBank("TimeBasedTrkg::TBSegments", segList.size());

        int[] hitIdxArray = new int[12];

        for (int i = 0; i < segList.size(); i++) {
            if (segList.get(i).get_Id() == -1) continue;

            for (int j = 0; j < hitIdxArray.length; j++) {
                hitIdxArray[j] = -1;
            }

            double chi2 = 0;

            bank.setShort("id",        i, (short) segList.get(i).get_Id());
            bank.setByte("superlayer", i, (byte) segList.get(i).get_Superlayer());
            bank.setByte("sector",     i, (byte) segList.get(i).get_Sector());
            if (TB) bank.setShort("status", i, (short) segList.get(i).get_Status());

            FittedCluster cls = segList.get(i).get_fittedCluster();
            bank.setShort("Cluster_ID",    i, (short) cls.get_Id());

            // bank.setFloat("avgWire",       i, (float) cls.getAvgwire());
            // bank.setByte ("size",          i, (byte)  segList.get(i).size());
            // bank.setFloat("fitSlope",      i, (float) cls.get_clusterLineFitSlope());
            // bank.setFloat("fitSlopeErr",   i, (float) cls.get_clusterLineFitSlopeErr());
            // bank.setFloat("fitInterc",     i, (float) cls.get_clusterLineFitIntercept());
            // bank.setFloat("fitIntercErr",  i, (float) cls.get_clusterLineFitInterceptErr());
            bank.setFloat("SegEndPoint1X", i, (float) segList.get(i).get_SegmentEndPoints()[0]);
            bank.setFloat("SegEndPoint1Z", i, (float) segList.get(i).get_SegmentEndPoints()[1]);
            bank.setFloat("SegEndPoint2X", i, (float) segList.get(i).get_SegmentEndPoints()[2]);
            bank.setFloat("SegEndPoint2Z", i, (float) segList.get(i).get_SegmentEndPoints()[3]);

            if (TB) {
                bank.setFloat("resiSum", i, (float) segList.get(i).get_ResiSum());
                bank.setFloat("timeSum", i, (float) segList.get(i).get_TimeSum());
            }

            // for (int j = 0; j < segList.get(i).size(); j++) {
            //     if (segList.get(i).get_Id() == -1) continue;
            //     if (j < hitIdxArray.length) {
            //         hitIdxArray[j] = segList.get(i).get(j).get_Id();
            //     }
            //     // Math.sqrt(12.) = 3.4641016151377544
            //     double residual = segList.get(i).get(j).get_ClusFitDoca() /
            //                       (segList.get(i).get(j).get_CellSize() / 3.4641016151377544);
            //     chi2 += residual * residual;
            // }
            // bank.setFloat("fitChisqProb", i,
            //               (float) ProbChi2perNDF.prob(chi2, segList.get(i).size() - 2));
            //
            // for (int j = 0; j < hitIdxArray.length; j++) {
            //     String hitStrg = "Hit";
            //     hitStrg += (j + 1);
            //     hitStrg += "_ID";
            //     bank.setShort(hitStrg, i, (short) hitIdxArray[j]);
            // }

            if (!TB) {
                bank.setByte("isOnTrack", i, (byte) (segList.get(i).isOnTrack ? 1 : 0));
                bank.setFloat("resiSum", i, (float) segList.get(i).get_ResiSum());
                bank.setFloat("timeSum", i, (float) segList.get(i).get_TimeSum());
                bank.setByte("status", i, (byte) segList.get(i).get_Status());
                bank.setShort("associatedCrossId", i, (short) segList.get(i).associatedCrossId);

                bank.setFloat("fitPlane_px", i, (float) segList.get(i).get_fitPlane().point().x());
                bank.setFloat("fitPlane_py", i, (float) segList.get(i).get_fitPlane().point().y());
                bank.setFloat("fitPlane_pz", i, (float) segList.get(i).get_fitPlane().point().z());
                bank.setFloat("fitPlane_nx", i, (float) segList.get(i).get_fitPlane().normal().x());
                bank.setFloat("fitPlane_ny", i, (float) segList.get(i).get_fitPlane().normal().y());
                bank.setFloat("fitPlane_nz", i, (float) segList.get(i).get_fitPlane().normal().z());
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
    private DataBank fillCrossesBank(DataEvent event, List<Cross> crossList, boolean TB) {
        if (TB && event.hasBank("TimeBasedTrkg::TBCrosses")) { // For second pass tracking
            ((HipoDataEvent) event).getHipoEvent().removeGroup("TimeBasedTrkg::TBCrosses");
        }
        if (!TB && event.hasBank("HitBasedTrkg::HBCrosses")) {
            ((HipoDataEvent) event).getHipoEvent().removeGroup("HitBasedTrkg::HBCrosses");
        }
        int banksize = 0;
        for (Cross cross : crossList) {
            if (cross.get_Id() != -1) banksize++;
        }

        DataBank bank;
        if (!TB) bank = event.createBank("HitBasedTrkg::HBCrosses",  banksize);
        else     bank = event.createBank("TimeBasedTrkg::TBCrosses", banksize);

        int idx = 0;
        for (Cross cross : crossList) {
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

            if (!TB) bank.setInt("recalc", idx, (int) cross.recalc);

            idx++;
        }

        return bank;
    }

    /**
     * Writes a list of tracks into the EvioEvent's bank.
     * @param event    the EvioEvent
     * @param candList the list of tracks
     * @param TB       boolean set to 1 if time-based and 0 otherwise.
     * @return         tracks bank
     */
    private DataBank fillTracksBank(DataEvent event, List<Track> candList, boolean TB,
                                    boolean DCHB2) {
        if (TB && event.hasBank("TimeBasedTrkg::TBTracks")) { // For second pass tracking
            ((HipoDataEvent) event).getHipoEvent().removeGroup("TimeBasedTrkg::TBTracks");
        }
        if (!TB && event.hasBank("HitBasedTrkg::HBTracks")) {
            ((HipoDataEvent) event).getHipoEvent().removeGroup("HitBasedTrkg::HBTracks");
        }

        DataBank bank;
        if (!TB) bank = event.createBank("HitBasedTrkg::HBTracks",  candList.size());
        else     bank = event.createBank("TimeBasedTrkg::TBTracks", candList.size());

        if (candList.size() == 0) return bank;

        for (int i = 0; i < candList.size(); i++) {
            bank.setShort("id",     i, (short) candList.get(i).get_Id());
            bank.setByte ("sector", i, (byte)  candList.get(i).get_Sector());
            bank.setByte ("q",      i, (byte)  candList.get(i).get_Q());
            if (!DCHB2) bank.setShort("status", i, (short) candList.get(i).get_Status());
            else        bank.setShort("status", i, (short) (100 + candList.get(i).get_Status() * 10
                                                            + candList.get(i).get_MissingSuperlayer()));

            if (candList.get(i).get_PreRegion1CrossPoint() != null) {
                bank.setFloat("c1_x",  i, (float) candList.get(i).get_PreRegion1CrossPoint().x());
                bank.setFloat("c1_y",  i, (float) candList.get(i).get_PreRegion1CrossPoint().y());
                bank.setFloat("c1_z",  i, (float) candList.get(i).get_PreRegion1CrossPoint().z());
                bank.setFloat("c1_ux", i, (float) candList.get(i).get_PreRegion1CrossDir().x());
                bank.setFloat("c1_uy", i, (float) candList.get(i).get_PreRegion1CrossDir().y());
                bank.setFloat("c1_uz", i, (float) candList.get(i).get_PreRegion1CrossDir().z());
            }
            if (candList.get(i).get_PostRegion3CrossPoint() != null) {
                bank.setFloat("c3_x",  i, (float) candList.get(i).get_PostRegion3CrossPoint().x());
                bank.setFloat("c3_y",  i, (float) candList.get(i).get_PostRegion3CrossPoint().y());
                bank.setFloat("c3_z",  i, (float) candList.get(i).get_PostRegion3CrossPoint().z());
                bank.setFloat("c3_ux", i, (float) candList.get(i).get_PostRegion3CrossDir().x());
                bank.setFloat("c3_uy", i, (float) candList.get(i).get_PostRegion3CrossDir().y());
                bank.setFloat("c3_uz", i, (float) candList.get(i).get_PostRegion3CrossDir().z());
            }
            if (candList.get(i).get_Region1TrackX() != null) {
                bank.setFloat("t1_x",  i, (float) candList.get(i).get_Region1TrackX().x());
                bank.setFloat("t1_y",  i, (float) candList.get(i).get_Region1TrackX().y());
                bank.setFloat("t1_z",  i, (float) candList.get(i).get_Region1TrackX().z());
                bank.setFloat("t1_px", i, (float) candList.get(i).get_Region1TrackP().x());
                bank.setFloat("t1_py", i, (float) candList.get(i).get_Region1TrackP().y());
                bank.setFloat("t1_pz", i, (float) candList.get(i).get_Region1TrackP().z());
            }

            bank.setFloat("pathlength", i, (float) candList.get(i).get_TotPathLen());
            bank.setFloat("Vtx0_x",     i, (float) candList.get(i).get_Vtx0().x());
            bank.setFloat("Vtx0_y",     i, (float) candList.get(i).get_Vtx0().y());
            bank.setFloat("Vtx0_z",     i, (float) candList.get(i).get_Vtx0().z());
            bank.setFloat("p0_x",       i, (float) candList.get(i).get_pAtOrig().x());
            bank.setFloat("p0_y",       i, (float) candList.get(i).get_pAtOrig().y());
            bank.setFloat("p0_z",       i, (float) candList.get(i).get_pAtOrig().z());
            bank.setFloat("chi2",       i, (float) candList.get(i).get_FitChi2());
            bank.setShort("ndf",        i, (short) candList.get(i).get_FitNDF());

            if (!TB) {
                bank.setFloat("_IntegralBdl",   i, (float) candList.get(i).get_IntegralBdl());
                bank.setFloat("_pathLength",    i, (float) candList.get(i).get_PathLength());
                bank.setFloat("_P",             i, (float) candList.get(i).get_P());
                bank.setByte ("fit_Successful", i, (byte) (candList.get(i).fit_Successful ? 1 : 0));
                bank.setShort("_missingSuperlayer", i, (short) candList.get(i).get_MissingSuperlayer());
                bank.setShort("_fitConvergenceStatus", i, (short) candList.get(i).get_FitConvergenceStatus());
                bank.setFloat("b_0", i, (float) candList.get(i).b[0]);
                bank.setFloat("b_1", i, (float) candList.get(i).b[1]);
                bank.setFloat("b_2", i, (float) candList.get(i).b[2]);

                bank.setShort("_StateVecAtReg1MiddlePlane", i,
                              (short) candList.get(i).get_StateVecAtReg1MiddlePlane().getId());

                if (candList.get(i).get_Trajectory() != null) {
                    bank.setShort("n_sv", i, (short) candList.get(i).get_Trajectory().size());
                    for (int it = 0; it < candList.get(i).get_Trajectory().size(); ++it)
                        bank.setShort("svid" + it, i, (short) candList.get(i).get_Trajectory().get(it).getId());
                }
            }

            if (!TB) {
                bank.setShort("Cross1_ID",  i, (short) candList.get(i).get(0).get_Id());
                bank.setShort("Cross2_ID",  i, (short) candList.get(i).get(1).get_Id());
                bank.setShort("Cross3_ID",  i, (short) candList.get(i).get(2).get_Id());
                continue;
            }

            // else
            if (candList.get(i).size() == 3) {
                bank.setShort("Cross1_ID", i, (short) candList.get(i).get(0).get_Id());
                bank.setShort("Cross2_ID", i, (short) candList.get(i).get(1).get_Id());
                bank.setShort("Cross3_ID", i, (short) candList.get(i).get(2).get_Id());
            }
            else if (candList.get(i).size() == 2) {
                bank.setShort("Cross1_ID", i, (short) candList.get(i).get(0).get_Id());
                bank.setShort("Cross2_ID", i, (short) candList.get(i).get(1).get_Id());
                bank.setShort("Cross3_ID", i, (short) -1);
            }
            else if (candList.get(i).size() == 1) {
                bank.setShort("Cross1_ID", i, (short) candList.get(i).get(0).get_Id());
                bank.setShort("Cross2_ID", i, (short) -1);
                bank.setShort("Cross3_ID", i, (short) -1);
            }
        }
        return bank;
    }

    /**
     * writes the covariance matrix from HB fits to be used for starting the Time-based tracking.
     * @param event    hipo event
     * @param candList list of tracks
     * @return         covariance matrix
     */
    private DataBank fillTrackCovMatBank(DataEvent event, List<Track> candList) {
        if (event.hasBank("TimeBasedTrkg::TBCovMat")) {
            ((HipoDataEvent) event).getHipoEvent().removeGroup("TimeBasedTrkg::TBCovMat");
        }

        if (event == null) return null;
        DataBank bank = event.createBank("TimeBasedTrkg::TBCovMat", candList.size());

        for (int i = 0; i < candList.size(); i++) {
            bank.setShort("id", i, (short) candList.get(i).get_Id());
            if (candList.get(i).get_CovMat() == null) continue;
            bank.setFloat("C11", i, (float) candList.get(i).get_CovMat().get(0, 0));
            bank.setFloat("C12", i, (float) candList.get(i).get_CovMat().get(0, 1));
            bank.setFloat("C13", i, (float) candList.get(i).get_CovMat().get(0, 2));
            bank.setFloat("C14", i, (float) candList.get(i).get_CovMat().get(0, 3));
            bank.setFloat("C15", i, (float) candList.get(i).get_CovMat().get(0, 4));
            bank.setFloat("C21", i, (float) candList.get(i).get_CovMat().get(1, 0));
            bank.setFloat("C22", i, (float) candList.get(i).get_CovMat().get(1, 1));
            bank.setFloat("C23", i, (float) candList.get(i).get_CovMat().get(1, 2));
            bank.setFloat("C24", i, (float) candList.get(i).get_CovMat().get(1, 3));
            bank.setFloat("C25", i, (float) candList.get(i).get_CovMat().get(1, 4));
            bank.setFloat("C31", i, (float) candList.get(i).get_CovMat().get(2, 0));
            bank.setFloat("C32", i, (float) candList.get(i).get_CovMat().get(2, 1));
            bank.setFloat("C33", i, (float) candList.get(i).get_CovMat().get(2, 2));
            bank.setFloat("C34", i, (float) candList.get(i).get_CovMat().get(2, 3));
            bank.setFloat("C35", i, (float) candList.get(i).get_CovMat().get(2, 4));
            bank.setFloat("C41", i, (float) candList.get(i).get_CovMat().get(3, 0));
            bank.setFloat("C42", i, (float) candList.get(i).get_CovMat().get(3, 1));
            bank.setFloat("C43", i, (float) candList.get(i).get_CovMat().get(3, 2));
            bank.setFloat("C44", i, (float) candList.get(i).get_CovMat().get(3, 3));
            bank.setFloat("C45", i, (float) candList.get(i).get_CovMat().get(3, 4));
            bank.setFloat("C51", i, (float) candList.get(i).get_CovMat().get(4, 0));
            bank.setFloat("C52", i, (float) candList.get(i).get_CovMat().get(4, 1));
            bank.setFloat("C53", i, (float) candList.get(i).get_CovMat().get(4, 2));
            bank.setFloat("C54", i, (float) candList.get(i).get_CovMat().get(4, 3));
            bank.setFloat("C55", i, (float) candList.get(i).get_CovMat().get(4, 4));
        }

        return bank;
    }

    /**
     * Writes a list of trajectories with data pulled from a list of tracks.
     * @param event  the EvioEvent
     * @param tracks the list of tracks
     * @return       trajectories bank
     */
    public DataBank fillTrajectoryBank(DataEvent event, List<Track> tracks) {
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

    private DataBank fillStateVecsBank(DataEvent event, List<Track> tracks) {
        int bankSize = 0;
        for (Track track : tracks) {
            bankSize += track.get_Trajectory().size();
            bankSize++;
        }

        if (event.hasBank("TimeBasedTrkg::StateVec")) {
            ((HipoDataEvent) event).getHipoEvent().removeGroup("TimeBasedTrkg::StateVec");
        }

        if (bankSize == 0) return null;
        DataBank bank = event.createBank("TimeBasedTrkg::StateVec", bankSize);

        List<StateVec> stateVecList = new ArrayList<>();
        for (Track track : tracks) {
            if (track == null)                  continue;
            if (track.get_Trajectory() == null) continue;

            for (int i = 0; i < track.get_Trajectory().size(); ++i) {
                stateVecList.add(track.get_Trajectory().get(i));
            }
            stateVecList.add(track.get_StateVecAtReg1MiddlePlane());
        }

        int idx = 0;

        for (StateVec stateVec : stateVecList) {
            bank.setShort("id",          idx, (short) stateVec.getId());
            bank.setFloat("x",           idx, (float) stateVec.x());
            bank.setFloat("y",           idx, (float) stateVec.y());
            bank.setFloat("thX",         idx, (float) stateVec.tanThetaX());
            bank.setFloat("thY",         idx, (float) stateVec.tanThetaY());
            bank.setFloat("z",           idx, (float) stateVec.getZ());
            bank.setFloat("b",           idx, (float) stateVec.getB());
            bank.setFloat("h",           idx, (float) stateVec.getProjector());
            bank.setFloat("_PathLength", idx, (float) stateVec.getPathLength());
            bank.setShort("_planeIdx",   idx, (short) stateVec.getPlaneIdx());

            idx++;
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
                             boolean             TB,
                             boolean             DCHB2) {

        if (event == null) return;
        if (fhits != null) event.appendBanks(rbw.fillHitsBank(event, fhits, TB));
        else return;

        if (clusters != null) event.appendBanks(rbw.fillClustersBank(event, clusters, TB, DCHB2));
        else return;

        if (segments != null) event.appendBanks(rbw.fillSegmentsBank(event, segments, TB));
        else return;

        if (crosses != null) event.appendBanks(rbw.fillCrossesBank(event, crosses, TB));
        else return;

        if (trkcands != null) {
            event.appendBanks(rbw.fillTracksBank(event, trkcands, TB, DCHB2));
            if (!TB) {
                event.appendBanks(rbw.fillTrackCovMatBank(event, trkcands));
                event.appendBanks(rbw.fillStateVecsBank(event, trkcands));
            }
            else {
                event.appendBanks(rbw.fillTrajectoryBank(event, trkcands));
            }
        }
        return;
    }

    public void fillAllHBBanks(DataEvent event, RecoBankWriter rbw,
                               List<FittedHit>     fhits,
                               List<FittedCluster> clusters,
                               List<Segment>       segments,
                               List<Cross>         crosses,
                               List<Track>         trkcands) {
        fillAllBanks(event, rbw, fhits, clusters, segments, crosses, trkcands, false, false);
        return;
    }

    public void fillAllHBBanksFinal(DataEvent event, RecoBankWriter rbw,
                                    List<FittedHit>     fhits,
                                    List<FittedCluster> clusters,
                                    List<Segment>       segments,
                                    List<Cross>         crosses,
                                    List<Track>         trkcands) {
        fillAllBanks(event, rbw, fhits, clusters, segments, crosses, trkcands, false, true);
    }

    public void fillAllTBBanks(DataEvent event, RecoBankWriter rbw,
                               List<FittedHit>     fhits,
                               List<FittedCluster> clusters,
                               List<Segment>       segments,
                               List<Cross>         crosses,
                               List<Track>         trkcands) {
        fillAllBanks(event, rbw, fhits, clusters, segments, crosses, trkcands, true, false);
        return;
    }
}
