package org.jlab.rec.dc.cross;

import java.util.ArrayList;
// TODO: Uncomment
// import org.jlab.clas.clas.math.FastMath;
import org.apache.commons.math3.util.FastMath;

import org.jlab.geom.prim.Point3D;
import org.jlab.detector.geant4.v2.DCGeant4Factory;

import org.jlab.rec.dc.hit.FittedHit;
import org.jlab.rec.dc.segment.Segment;

/**
 * The crosses are objects used to find tracks and are characterized by a 3-D
 * point and a direction unit vector.
 *
 * @author ziegler
 */
public class Cross extends ArrayList<Segment> implements Comparable<Cross> {

    private static final long serialVersionUID = 5317526429163382618L;

    private int _Sector; // sector[1...6]
    private int _Region; // region[1,...3]
    private int _Id;	 // cross Id

    // Point parameters:
    private Point3D _Point;    // A 3-D point characterizing the position of the cross in the tilted
                               //     coordinate system.
    private Point3D _PointErr; // A 3-dimensional error on the 3-D point characterizing the position
                               //     of the cross in the tilted coordinate system.
    private Point3D _Dir;      // The cross unit direction vector.
    private Point3D _DirErr;   // The cross unit direction vector's error.

    private Segment _seg1; // The segment in the first superlayer, used to make a cross.
    private Segment _seg2; // The segment in the second superlayer, used to make a cross.
    public boolean isPseudoCross = false;

    // Math.sin(Math.toRadians(25.)) = 0.42261826174069944
    // Math.cos(Math.toRadians(25.)) = 0.9063077870366499
    private double cos_tilt = 0.9063077870366499;
    private double sin_tilt = 0.42261826174069944;

    public int recalc;

    /**
     * Constructor.
     * @param sector the sector (1...6)
     * @param region the region (1...3)
     * @param rid    the cross ID (if there are only 3 crosses in the event, the ID corresponds to
     *               the region index
     */
    public Cross(int sector, int region, int rid) {
        this._Sector = sector;
        this._Region = region;
        this._Id = rid;
    }

    public static long getSerialversionuid() {return serialVersionUID;}

    public int get_Sector() {return _Sector;}
    public void set_Sector(int _Sector) {this._Sector = _Sector;}

    public int get_Region() {return _Region;}
    public void set_Region(int _Region) {this._Region = _Region;}

    public int get_Id() {return _Id;}
    public void set_Id(int _Id) {this._Id = _Id;}

    public Point3D get_Point() {return _Point;}
    public void set_Point(Point3D _Point) {this._Point = _Point;}

    public Point3D get_PointErr() {return _PointErr;}
    public void set_PointErr(Point3D _PointErr) {this._PointErr = _PointErr;}

    public Point3D get_Dir() {return _Dir;}
    public void set_Dir(Point3D _Dir) {this._Dir = _Dir;}

    public Point3D get_DirErr() {return _DirErr;}
    public void set_DirErr(Point3D _DirErr) {this._DirErr = _DirErr;}

    public Segment get_Segment1() {return _seg1;}
    public void set_Segment1(Segment seg1) {this._seg1 = seg1;}

    public Segment get_Segment2() {return _seg2;}
    public void set_Segment2(Segment seg2) {this._seg2 = seg2;}

    /**
     * Compares crosses by azimuth angle values.
     * @param arg the cross to which this one is compared
     * @return    integer describing which int has the higher azimuth angle value
     */
    @Override
    public int compareTo(Cross arg) {
        double theta_rad1 = FastMath.atan2(this.get_Point().y(), this.get_Point().x());
        double theta_deg1 = (theta_rad1 / Math.PI * 180.) + (theta_rad1 > 0 ? 0 : 360.);
        double theta_rad2 = FastMath.atan2(arg.get_Point().y(), arg.get_Point().x());
        double theta_deg2 = (theta_rad2 / Math.PI * 180.) + (theta_rad2 > 0 ? 0 : 360.);

        if (theta_deg1 < theta_deg2) return 1;
        else                         return -1;
    }

    /**
     * Sets the position and direction unit vector for the cross.
     * @param DcDetector DC detector geometry
     */
    public void set_CrossParams(DCGeant4Factory DcDetector) {

        // System.out.println("[set_CrossParams] 00: init");
        double z = DcDetector.getRegionMidpoint(this.get_Region() - 1).z;

        // Math.cos(Math.toRadians(6.)) / Math.sin(Math.toRadians(6.)) = 9.514364454222585
        double wy_over_wx = 9.514364454222585;
        double val_sl1 = this._seg1.get_fittedCluster().get_clusterLineFitSlope();
        double val_sl2 = this._seg2.get_fittedCluster().get_clusterLineFitSlope();
        double val_it1 = this._seg1.get_fittedCluster().get_clusterLineFitIntercept();
        double val_it2 = this._seg2.get_fittedCluster().get_clusterLineFitIntercept();

        double x = 0.5 * (val_it1 + val_it2) + 0.5 * z * (val_sl1 + val_sl2);
        double y = 0.5 * wy_over_wx * (val_it2 - val_it1) + 0.5 * wy_over_wx * z * (val_sl2 - val_sl1);

        this.set_Point(new Point3D(x, y, z));

        double tanThX = val_sl2;
        double tanThY = FastMath.atan2(y, z);
        double uz = 1. / Math.sqrt(1 + tanThX * tanThX + tanThY * tanThY);
        double ux = uz * tanThX;
        double uy = uz * tanThY;

        this.set_Dir(new Point3D(ux, uy, uz));

        if (this.get_Dir().z() == 0) {
            // System.out.println("[set_CrossParams] 01A: cross' dir.z == 0.");
            return;
        }
        // Error ...  propagated from errors on slopes and intercepts
        double err_sl1  = this._seg1.get_fittedCluster().get_clusterLineFitSlopeErr();
        double err_sl2  = this._seg2.get_fittedCluster().get_clusterLineFitSlopeErr();
        double err_it1  = this._seg1.get_fittedCluster().get_clusterLineFitInterceptErr();
        double err_it2  = this._seg2.get_fittedCluster().get_clusterLineFitInterceptErr();
        double err_cov1 = this._seg1.get_fittedCluster().get_clusterLineFitSlIntCov();
        double err_cov2 = this._seg2.get_fittedCluster().get_clusterLineFitSlIntCov();

        double err_x_fix = 0.5 * Math.sqrt(err_it1 * err_it1 + err_it2 * err_it2
                                           + z * z * (err_sl1 * err_sl1 + err_sl2 * err_sl2)
                                           + 2 * z * err_cov1 + 2 * z * err_cov2);
        double err_y_fix = 0.5 * wy_over_wx * Math.sqrt(err_it1 * err_it1 + err_it2 * err_it2
                                           + z * z * (err_sl1 * err_sl1 + err_sl2 * err_sl2)
                                           + 2 * z * err_cov1 + 2 * z * err_cov2);

        this.set_PointErr(new Point3D(err_x_fix, err_y_fix, 0));

        double inv_N_sq = 1. / (0.25 * val_sl1 * val_sl1 * (1 + wy_over_wx * wy_over_wx)
                        + 0.25 * val_sl2 * val_sl2 * (1 + wy_over_wx * wy_over_wx)
                        + 0.5  * val_sl1 * val_sl2 * (1 - wy_over_wx * wy_over_wx) + 1);

        double inv_N = Math.sqrt(inv_N_sq);
        double del_inv_n_del_sl1err = -0.25 * (val_sl1 * (1 + wy_over_wx * wy_over_wx)
                                             + val_sl2 * (1 - wy_over_wx * wy_over_wx))
                                                       * inv_N * inv_N_sq;

        double del_inv_n_del_sl2err = -0.25 * (val_sl2 * (1 + wy_over_wx * wy_over_wx)
                                             + val_sl1 * (1 - wy_over_wx * wy_over_wx))
                                                       * inv_N * inv_N_sq;

        double err_x2 = 0.5 * (inv_N + (val_sl1 + val_sl2) * del_inv_n_del_sl1err)
                            * (inv_N + (val_sl1 + val_sl2) * del_inv_n_del_sl1err)
                            * err_sl1 * err_sl1
                      + 0.5 * (inv_N + (val_sl1 + val_sl2) * del_inv_n_del_sl2err)
                            * (inv_N + (val_sl1 + val_sl2) * del_inv_n_del_sl2err)
                            * err_sl2 * err_sl2;

        double err_y2 = 0.5 * (-inv_N * wy_over_wx
                               + (val_sl2 - val_sl1) * wy_over_wx * del_inv_n_del_sl1err)
                            * (-inv_N * wy_over_wx
                               + (val_sl2 - val_sl1) * wy_over_wx * del_inv_n_del_sl1err)
                            * err_sl1 * err_sl1
                      + 0.5 * (inv_N * wy_over_wx + (val_sl2 - val_sl1) * del_inv_n_del_sl2err)
                            * (inv_N * wy_over_wx + (val_sl2 - val_sl1) * del_inv_n_del_sl2err)
                            * err_sl2 * err_sl2;

        double err_z2 = del_inv_n_del_sl1err * del_inv_n_del_sl1err * err_sl1 * err_sl1
                      + del_inv_n_del_sl2err * del_inv_n_del_sl2err * err_sl2 * err_sl2;

        double err_xDir = Math.sqrt(err_x2);
        double err_yDir = Math.sqrt(err_y2);
        double err_zDir = Math.sqrt(err_z2);

        this.set_DirErr(new Point3D(err_xDir, err_yDir, err_zDir));

        if (this._seg1.get_Id() == -1 || this._seg2.get_Id() == -1) this.isPseudoCross = true;

        // System.out.println("[set_CrossParams] 01B cross parameters recalculated successfully.");
    }

    /**
     * Gets rotated coordinates from tilted sector coordinate system to the sector's coordinate
     * system.
     * @param X NOTE: Parameter lacks description
     * @param Y NOTE: Parameter lacks description
     * @param Z NOTE: Parameter lacks description
     * @return  point3D containing the rotated coordinates
     */
    public Point3D getCoordsInSector(double X, double Y, double Z) {
        return new Point3D(X * cos_tilt + Z * sin_tilt, Y, -X * sin_tilt + Z * cos_tilt);
    }

    /**
     * Gets rotates coordinates from tilted sector coordinate system to the lab frame.
     * @param X NOTE: Parameter lacks description
     * @param Y NOTE: Parameter lacks description
     * @param Z NOTE: Parameter lacks description
     * @return  Point3D containing the rotated coordinates
     */
    public Point3D getCoordsInLab(double X, double Y, double Z) {
        Point3D PointInSec = this.getCoordsInSector(X, Y, Z);
        // Math.toRadians(60.) = 1.0471975511965976
        double rx = PointInSec.x() * FastMath.cos((this.get_Sector() - 1) * 1.0471975511965976)
                  - PointInSec.y() * FastMath.sin((this.get_Sector() - 1) * 1.0471975511965976);
        double ry = PointInSec.x() * FastMath.sin((this.get_Sector() - 1) * 1.0471975511965976)
                  + PointInSec.y() * FastMath.cos((this.get_Sector() - 1) * 1.0471975511965976);

        return new Point3D(rx, ry, PointInSec.z());
    }

    // NOTE: Method lacks JavaDoc description
    public Point3D getCoordsInTiltedSector(double X, double Y, double Z) {
        // Math.toRadians(-60.) = -1.0471975511965976
        double rx = X * FastMath.cos((this.get_Sector() - 1) * -1.0471975511965976)
                  - Y * FastMath.sin((this.get_Sector() - 1) * -1.0471975511965976);
        double ry = X * FastMath.sin((this.get_Sector() - 1) * -1.0471975511965976)
                  + Y * FastMath.cos((this.get_Sector() - 1) * -1.0471975511965976);

        return new Point3D(rx * cos_tilt - Z * sin_tilt, ry, rx * sin_tilt + Z * cos_tilt);
    }

    // NOTE: Method lacks JavaDoc description
    public void set_CrossDirIntersSegWires() {
        double wy_over_wx = (FastMath.cos(Math.toRadians(6.)) / FastMath.sin(Math.toRadians(6.)));
        double val_sl1 = this._seg1.get_fittedCluster().get_clusterLineFitSlope();
        double val_sl2 = this._seg2.get_fittedCluster().get_clusterLineFitSlope();
        double val_it1 = this._seg1.get_fittedCluster().get_clusterLineFitIntercept();
        double val_it2 = this._seg2.get_fittedCluster().get_clusterLineFitIntercept();

        for(int i = 0; i < this.get_Segment1().size(); i++) {
            this.calc_IntersectPlaneAtZ(this.get_Segment1().get(i).get_Z(), wy_over_wx,
                                        val_sl1, val_sl2, val_it1, val_it2,
                                        this.get_Segment1().get(i));
        }
        for(int i =0; i<this.get_Segment2().size(); i++) {
            this.calc_IntersectPlaneAtZ(this.get_Segment2().get(i).get_Z(), wy_over_wx,
                                        val_sl1, val_sl2, val_it1, val_it2,
                                        this.get_Segment2().get(i));
        }
    }

    // NOTE: Method lacks JavaDoc description
    private void calc_IntersectPlaneAtZ(double z,       double wy_over_wx,
                                        double val_sl1, double val_sl2,
                                        double val_it1, double val_it2,
                                        FittedHit hit) {

        double x = 0.5 * (val_it1 + val_it2) + 0.5 * z * (val_sl1 + val_sl2);
        double y = 0.5 * wy_over_wx * (val_it2 - val_it1)
                 + 0.5 * wy_over_wx * z * (val_sl2 - val_sl1);

        hit.setCrossDirIntersWire(new Point3D(x,y,z));
    }

    /**
     * Writes and returns a string containing the cross' information.
     * @return cross' information encoded in a string
     */
    public String getInfo() {
        return "DC cross:" +
               "\nID     : " + this.get_Id() +
               "\nSector : " + this.get_Sector() +
               "\nRegion : " + this.get_Region() +
               "\nPoint  : " + this.get_Point().toString() +
               "\nDir    : " + this.get_Dir().toString();
    }

    /**
     * Writes and returns a string containing all of the cross' data.
     * @return cross' data encoded in a string
     */
    public String getDetailedInfo() {
        return "DC Cross " + this.get_Id() + ":" +
               "\n  Sector         : " + this.get_Sector() +
               "\n  Region         : " + this.get_Region() +
               "\n  Point          : " + this.get_Point() +
               "\n  Point Err      : " + this.get_PointErr() +
               "\n  Direction      : " + this.get_Dir() +
               "\n  Direcction Err : " + this.get_DirErr() +
               "\n--------------------------------------------\n";
    }
}
