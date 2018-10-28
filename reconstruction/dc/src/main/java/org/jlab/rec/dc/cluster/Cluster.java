package org.jlab.rec.dc.cluster;

import java.util.ArrayList;

import org.jlab.rec.dc.hit.Hit;

/**
 * A cluster in the DC consists of an array of hits that are grouped together according to the
 * clustering algorithm ubicated in the ClusterFinder class.
 *
 * @author ziegler
 */
public class Cluster extends ArrayList<Hit> {

    private static final long serialVersionUID = 9153980362683755204L;

    private int _Sector;     //	sector [1...6]
    private int _Superlayer; //	superlayer [1,...6]
    private int _Id;		 //	cluster ID

    /**
     * Constructor method.
     * @param sector     the sector (1...6)
     * @param superlayer the superlayer (1...6)
     * @param cid        the cluster ID, an incremental integer corresponding to the cluster formed
     *                   in the series of clusters
     */
    public Cluster(int sector, int superlayer, int cid) {
        this._Sector = sector;
        this._Superlayer = superlayer;
        this._Id = cid;
    }

    /**
     * Creates a new list of hits characterized by its sector, superlayer and ID number.
     * @param hit the first hit in the list of hits composing the cluster
     * @param cid the id of the cluster
     * @return    the cluster
     */
    public Cluster newCluster(Hit hit, int cid) {
        return new Cluster(hit.get_Sector(), hit.get_Superlayer(), cid);
    }

    public int get_Sector() {return _Sector;}
    public void set_Sector(int _Sector) {this._Sector = _Sector;}

    public int get_Superlayer() {return _Superlayer;}
    public void set_Superlayer(int _Superlayer) {this._Superlayer = _Superlayer;}

    /** @return region (1...3) */
    public int get_Region() {return (int) (this._Superlayer + 1) / 2;}
    /** @return superlayer 1 or 2 in region (1...3) */
    public int get_RegionSlayer() {return (this._Superlayer + 1) % 2 + 1;}

    public int get_Id() {return _Id;}
    public void set_Id(int _Id) {this._Id = _Id;}

    /**
     * Writes and returns a string containing the cluster's information, including the locations
     * and the number of hits contained in it.
     * @return cluster's information encoded in a string
     */
    public String printInfo() {
        return "DC cluster:" +
               "\nID         : " + this.get_Id() +
               "\nSector     : " + this.get_Sector() +
               "\nSuperlayer : " + this.get_Superlayer() +
               "\nSize       : " + this.size();
    }
}
