package org.jlab.rec.dc.cluster;

import java.util.ArrayList;

import org.jlab.geom.prim.Line3D;
import org.jlab.rec.dc.hit.FittedHit;

/**
 * A fitted cluster in the DC consists of an array of hits that are grouped together according to
 * the algorithm of the ClusterFinder class and have been fit using the wire position information
 * and subsequent time-based information at the midplane.
 *
 * @author ziegler
 */
public class FittedCluster extends ArrayList<FittedHit> implements Comparable<FittedCluster> {

    private static final long serialVersionUID = 7240609802152999866L;

    private int _Sector;     // sector[1...6]
    private int _Superlayer; // superlayer [1,...6]
    private int _Id;		 //	cluster Id, which is the index in the sequence of formed clusters.

    private Line3D _clusLine;    // The line corresponding to the linear fit to the cluster's hits.
                                 //     The line is defined as a point on the line and a unit
                                 //     direction vector along the line.
    private Line3D _clusLineErr; // The cluster line error.

    private double _fitProb = -1; // the linear fit chi^2 probability
    private double _Chisq = Double.POSITIVE_INFINITY;

    private double _clusterLineFitSlope;        // The slope of the line fitted to the cluster.
    private double _clusterLineFitintercept;    // The intercept of the line fitted to the cluster.
    private double _clusterLineFitSlopeErr;     // The error in the slope.
    private double _clusterLineFitinterceptErr; // The error in the intercept.
    private double _clusterLineFitSlIntCov;

    private double _clusterLineFitSlopeMP;
    private double _clusterLineFitInterceptMP;
    private double _clusterLineFitSlopeErrMP;
    private double _clusterLineFitInterceptErrMP;

    private int[][] _Status;

    /**
     * Constructs a FittedCluster using the data from a cluster.
     * @param rawCluster a Cluster fit using hit-based tracking information
     */
    public FittedCluster(Cluster rawCluster) {
        this._Sector     = rawCluster.get_Sector();
        this._Superlayer = rawCluster.get_Superlayer();
        this._Id         = rawCluster.get_Id();

        // Add the hits to the defined cluster
        for (int i = 0; i < rawCluster.size(); i++) {
            FittedHit fhit = new FittedHit(rawCluster.get(i).get_Sector(),
                                           rawCluster.get(i).get_Superlayer(),
                                           rawCluster.get(i).get_Layer(),
                                           rawCluster.get(i).get_Wire(),
                                           rawCluster.get(i).get_TDC(),
                                           rawCluster.get(i).get_Id());
            fhit.set_DocaErr(rawCluster.get(i).get_DocaErr());
            fhit.set_CellSize(rawCluster.get(i).get_CellSize());
            fhit.set_Id(rawCluster.get(i).get_Id());
            this.add(fhit);
        }
    }

    /** Explicit constructor for a FittedCluster. */
    public FittedCluster(int sector, int superlayer, int id, ArrayList<FittedHit> fHits) {
        this._Sector = sector;
        this._Superlayer = superlayer;
        this._Id = id;

        for (int i = 0; i < fHits.size(); i++) this.add(fHits.get(i));
    }


    public int get_Sector() {return _Sector;}
    public void set_Sector(int _Sector) {this._Sector = _Sector;}

    public int get_Superlayer() {return _Superlayer;}
    public void set_Superlayer(int _Superlayer) {this._Superlayer = _Superlayer;}

    /** @return region (1...3) */
    public int get_Region() {return (int) (this._Superlayer + 1) / 2;}

    /** @return superlayer 1 or 2 in region */
    public int get_RegionSlayer() {return (this._Superlayer + 1) % 2 + 1;}

    public int get_Id() {return _Id;}
    public void set_Id(int _Id) {this._Id = _Id;}

    public int[][] get_Status() {return _Status;}
    public void set_Status(int[][] _Status) {this._Status = _Status;}

    public Line3D get_clusLine() {return _clusLine;}
    public void set_clusLine(Line3D _clusLine) {this._clusLine = _clusLine;}

    public Line3D get_clusLineErr() {return _clusLineErr;}
    public void set_clusLineErr(Line3D _clusLineErr) {this._clusLineErr = _clusLineErr;}

    public double get_fitProb() {return _fitProb;}
    public void set_fitProb(double _fitChisq) {this._fitProb = _fitChisq;}

    public double get_Chisq() {return _Chisq;}
    public void set_Chisq(double _Chisq) {this._Chisq = _Chisq;}

    public double get_clusterLineFitSlope() {return _clusterLineFitSlope;}
    public void set_clusterLineFitSlope(double _clusterLineFitSlope) {
        this._clusterLineFitSlope = _clusterLineFitSlope;
    }

    public double get_clusterLineFitSlopeErr() {return _clusterLineFitSlopeErr;}
    public void set_clusterLineFitSlopeErr(double _clusterLineFitSlopeErr) {
        this._clusterLineFitSlopeErr = _clusterLineFitSlopeErr;
    }

    public double get_clusterLineFitIntercept() {return _clusterLineFitintercept;}
    public void set_clusterLineFitIntercept(double _clusterLineFitIntercept) {
        this._clusterLineFitintercept = _clusterLineFitIntercept;
    }

    public double get_clusterLineFitInterceptErr() {return _clusterLineFitinterceptErr;}
    public void set_clusterLineFitInterceptErr(double _clusterLineFitinterceptErr) {
        this._clusterLineFitinterceptErr = _clusterLineFitinterceptErr;
    }

    public double get_clusterLineFitSlIntCov() {return _clusterLineFitSlIntCov;}
    public void set_clusterLineFitSlIntCov(double _clusterLineFitSlIntCov) {
        this._clusterLineFitSlIntCov = _clusterLineFitSlIntCov;
    }

    public double get_clusterLineFitSlopeMP() {return _clusterLineFitSlopeMP;}
    public void set_clusterLineFitSlopeMP(double _clusterLineFitSlopeMP) {
        this._clusterLineFitSlopeMP = _clusterLineFitSlopeMP;
    }

    public double get_clusterLineFitSlopeErrMP() {return _clusterLineFitSlopeErrMP;}
    public void set_clusterLineFitSlopeErrMP(double _clusterLineFitSlopeErrMP) {
        this._clusterLineFitSlopeErrMP = _clusterLineFitSlopeErrMP;
    }

    public double get_clusterLineFitInterceptMP() {return _clusterLineFitInterceptMP;}
    public void set_clusterLineFitInterceptMP(double _clusterLineFitInterceptMP) {
        this._clusterLineFitInterceptMP = _clusterLineFitInterceptMP;
    }

    public double get_clusterLineFitInterceptErrMP() {return _clusterLineFitInterceptErrMP;}
    public void set_clusterLineFitInterceptErrMP(double _clusterLineFitInterceptErrMP) {
        this._clusterLineFitInterceptErrMP = _clusterLineFitInterceptErrMP;
    }

    /** @return Average wire position in a cluster. */
    public double getAvgwire() {
        double avewire = 0;
        int hSize = this.size();
        for (int h = 0; h < hSize; h++) avewire += this.get(h).get_Wire();
        return ((double) avewire / hSize);
    }

    /** @return the cluster from which the fitted cluster was created */
    public Cluster getBaseCluster() {
        return new Cluster(this.get_Sector(), this.get_Superlayer(), this.get_Id());
    }

    /**
     * Cluster comparator based on number of hits in the cluster.
     * @return an int describing which cluster is bigger.
     */
    @Override
    public int compareTo(FittedCluster o) {
        if (this.size() > o.size()) return 1;
        else return 0;
    }

    /**
     * Gets the current tracking status for the fitted cluster.
     * @return the tracking status: -1: no fits yet
     *                               0: hit-based tracking done
     *                               1: time-based tracking done
     */
    public int get_TrkgStatus() {

        boolean postHitBasedFit  = true;
        boolean postTimeBasedFit = true;

        for (FittedHit fhit : this) {
            if (fhit.get_TrkgStatus() <= 0)  postTimeBasedFit = false;
            if (fhit.get_TrkgStatus() == -1) postHitBasedFit = false;
        }
        if      (postTimeBasedFit) return 1;
        else if (postHitBasedFit)  return 0;
        else                       return -1;
    }

    /**
     * Writes and returns a string containing the fitted cluster's information.
     * @return fitted cluster's information encoded in a string
     */
    public String printInfo() {
        return "Fitted DC cluster:" +
               "\nID         : " + this.get_Id() +
               "\nSector     : " + this.get_Sector() +
               "\nSuperlayer : " + this.get_Superlayer() +
               "\nSize       : " + this.size() +
               "\nfit chi2   : " + this.get_fitProb();
    }
}
