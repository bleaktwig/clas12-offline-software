package org.jlab.rec.dc.segment;

import java.util.ArrayList;
import java.util.List;

import org.jlab.geom.prim.Plane3D;
import org.jlab.geom.prim.Point3D;
import org.jlab.geom.prim.Vector3D;
import org.jlab.detector.geant4.v2.DCGeant4Factory;

import org.jlab.rec.dc.Constants;
import org.jlab.rec.dc.hit.FittedHit;
import org.jlab.rec.dc.cluster.FittedCluster;
import org.jlab.rec.dc.trajectory.SegmentTrajectory;

/**
 * A class to describe segment objects, where a Segment is a fitted cluster that has been pruned of
 * hits with bad residuals
 *
 * @author ziegler
 */
public class Segment extends ArrayList<FittedHit> implements Comparable<Segment> {

    private static final long serialVersionUID = -997960312423538455L;
    private FittedCluster _fittedCluster;
    public boolean isOnTrack = false;

    private int     _Sector;     // sector [1...6]
    private int     _Superlayer; // superlayer [1,...6]
    private int     _Id;		 // associated cluster Id, corresponds to the index in the sequence
                                 //     of found segments.
    private double  _ResiSum;    // sum of residuals for hits in segment
    private double  _TimeSum;    // sum of corrected (T0-subtracted) times for all hits in segment
    private Plane3D _fitPlane;   // the plane containing the segment fitted-line representation

    private SegmentTrajectory _Trajectory;
    private double[] _SegmentEndPoints;

    private int _Status = 1;     // not yet implemented
    public int associatedCrossId = -1;

    /**
     * Construct the segment from a fitted cluster.
     * @param fCluster the fitted Cluster
     */
    public Segment(FittedCluster fCluster) {
        for (int i = 0; i < fCluster.size(); i++) this.add(fCluster.get(i));
        this.set_fittedCluster(fCluster);
        this._Sector     = fCluster.get_Sector();
        this._Superlayer = fCluster.get_Superlayer();
        this._Id         = fCluster.get_Id();
        this.set_Status(Status());
    }

    /** NOTE: Lacks description */
    public int Status() {
        int stat = 0;

        int L[] = new int[6];
        for (int l = 0; l < this.size(); l++) L[this.get(l).get_Layer() - 1]++;
        for (int l = 0; l < 6; l++) {
            if (L[l] == 0 || L[l] > 2) stat = 1;
        }
        return stat;
    }

    public FittedCluster get_fittedCluster() {return _fittedCluster;}
    public void set_fittedCluster(FittedCluster _fittedCluster) {
        this._fittedCluster = _fittedCluster;
    }

    public int get_Sector() {return _Sector;}
    public void set_Sector(int _Sector) {this._Sector = _Sector;}

    public int get_Superlayer() {return _Superlayer;}
    public void set_Superlayer(int _Superlayer) {this._Superlayer = _Superlayer;}

    public int get_Id() {return _Id;}
    public void set_Id(int _Id) {this._Id = _Id;}

    public int get_Region() {return (this._Superlayer + 1) / 2;}
    public int get_RegionSlayer() {return (this._Superlayer + 1) % 2 + 1;}

    public double get_ResiSum() {return _ResiSum;}
    public void set_ResiSum(double _ResiSum) {this._ResiSum = _ResiSum;}

    public double get_TimeSum() {return _TimeSum;}
    public void set_TimeSum(double _TimeSum) {this._TimeSum = _TimeSum;}

    public Plane3D get_fitPlane() {return _fitPlane;}

    public SegmentTrajectory get_Trajectory() {return _Trajectory;}
    public void set_Trajectory(SegmentTrajectory _Trajectory) {this._Trajectory = _Trajectory;}

    public double[] get_SegmentEndPoints() {return _SegmentEndPoints;}
    public void set_SegmentEndPoints(double[] _SegmentEndPoints) {
        this._SegmentEndPoints = _SegmentEndPoints;
    }

    public int get_Status() {return _Status;}
    public void set_Status(int _Status) {this._Status = _Status;}

    /**
     * Checks if another segment is close to this one.
     * @param otherseg matching cluster in other superlayer in a region
     * @return         a region-segment proximity condition
     */
    public boolean isCloseTo(Segment otherseg) {
        /*
        A region-segment contains two segments if they are in the same sector and region and satisfy
        the proximity condition:
        |Xwires_2 - Xwires_1| = a * Xwires1 + b
        where a and b are DC parameters set by DC_RSEG_a and DC_RSEG_B.
        */

        boolean value = false;

        if (Math.abs(this.getAvgwire() - otherseg.getAvgwire()) <
                Constants.DC_RSEG_A * this.getAvgwire() + Constants.DC_RSEG_B) {
            value = true;
        }

        return value;
    }

    // NOTE: Lacks JavaDoc comment
    public boolean hasNoMatchingSegment(List<Segment> othersegs) {
        /*
        A region-segment contains two segments if they are in the same sector and region and satisfy
        the proximity condition:
        |Xwires_2 - Xwires_1| = a * Xwires_1 + b
        where a and b are DC parameters set by DC_RSEG_a and DC_RSEG_B.
        */
        boolean value = true;

        for(Segment otherseg : othersegs) {
            if (Math.abs(this.getAvgwire() - otherseg.getAvgwire()) <
                    Constants.DC_RSEG_A * this.getAvgwire() + Constants.DC_RSEG_B) {
                // Found a matching segment
                value = false;
            }
        }

        return value;
    }

    // NOTE: Lacks JavaDoc comment
    public boolean hasConsistentSlope(Segment otherseg) {
        if (this.get_fitPlane() != null && otherseg.get_fitPlane() != null &&
                Math.abs(Math.toDegrees(Math.acos(this.get_fitPlane().normal().dot(otherseg.get_fitPlane().normal()))) - 12.) < Constants.SEGMENTPLANESANGLE) {
            return true;
        }
        return false;
    }

    /**
     * Gets the average wire number for the segment, used in the proximity condition employed in
     * segment matching.
     * @return the average wire number for the segment
     */
    public double getAvgwire() {
        double avewire = 0;
        int hSize = this.size();
        for (int h = 0; h < hSize; h++) avewire += this.get(h).get_Wire();

        return ((double) avewire / hSize);
    }

    /**
     * Sets the segment endpoints in the sector coordinate system for ced display.
     * @param DcDetector DC detector geometry
     */
    public void set_SegmentEndPointsSecCoordSys(DCGeant4Factory DcDetector) {

        double Z_1 = DcDetector.getWireMidpoint(this.get_Sector()-1, this.get_Superlayer()-1, 0, 0).z;
        double X_1 = this.get_fittedCluster().get_clusterLineFitSlope() * Z_1 +
                     this.get_fittedCluster().get_clusterLineFitIntercept();

        // Math.sin(Math.toRadians(25.)) = 0.42261826174069944
        // Math.cos(Math.toRadians(25.)) = 0.9063077870366499
        double x1 =  0.9063077870366499  * X_1 + 0.42261826174069944 * Z_1;
        double z1 = -0.42261826174069944 * X_1 + 0.9063077870366499  * Z_1;

        double Z_2 = DcDetector.getWireMidpoint(this.get_Sector()-1, this.get_Superlayer()-1, 5, 0).z;
        double X_2 = this.get_fittedCluster().get_clusterLineFitSlope() * Z_2 +
                     this.get_fittedCluster().get_clusterLineFitIntercept();

        double x2 =  0.9063077870366499  * X_2 + 0.42261826174069944 * Z_2;
        double z2 = -0.42261826174069944 * X_2 + 0.9063077870366499  * Z_2;

        double[] EndPointsArray = new double[4];
        EndPointsArray[0] = x1;
        EndPointsArray[1] = z1;
        EndPointsArray[2] = x2;
        EndPointsArray[3] = z2;

        this.set_SegmentEndPoints(EndPointsArray);
    }

    /**
     * Sets the plane containing the segment fitted-line representation.
     * @param DcDetector DC detector geometry
     */
    public void set_fitPlane(DCGeant4Factory DcDetector) {
        if (this.get_fittedCluster().get_clusLine() == null) return;

        this.set_SegmentEndPointsSecCoordSys(DcDetector);
        double dir_x = this.get_fittedCluster().get_clusLine().end().x()
                     - this.get_fittedCluster().get_clusLine().origin().x();
        double dir_y = this.get_fittedCluster().get_clusLine().end().y()
                     - this.get_fittedCluster().get_clusLine().origin().y();
        double dir_z = this.get_fittedCluster().get_clusLine().end().z()
                     - this.get_fittedCluster().get_clusLine().origin().z();

        double dir = Math.sqrt(dir_x * dir_x + dir_y * dir_y + dir_z * dir_z);

        dir_x /= dir;
        dir_y /= dir;
        dir_z /= dir;

        Vector3D dirVec = new Vector3D(dir_x, dir_y, dir_z);

        this._fitPlane = calc_fitPlane(this.get_fittedCluster().get_clusLine().origin(), dirVec);
    }

    /**
     * Calculates a plane containing the segment line, characterized by a point on the plane and a
     * unit normal vector.
     * @param refPoint the reference point on the segment plane
     * @param refDir   the normal vector on the segment plane
     * @return         the plane containing the segment line
     */
    private Plane3D calc_fitPlane(Point3D refPoint, Vector3D refDir) {

        // Math.sin(Math.toRadians(6)) = 0.10452846326765346
        // Math.cos(Math.toRadians(6)) = 0.9945218953682733
        double X = Math.pow(-1, (this.get_Superlayer() - 1)) * 0.10452846326765346;
        double Y = 0.9945218953682733;

        Vector3D plDir = new Vector3D(X, Y, 0);
        Vector3D normDir = plDir.cross(refDir);

        if (normDir.mag() > 1.e-10) normDir.scale(1. / normDir.mag());
        else                        return new Plane3D(new Point3D(0, 0, 0), new Vector3D(0, 0, 0));

        return new Plane3D(refPoint, normDir);
    }

    @Override
    /**
     * Compares two segments.
     * @param arg segment to which this one is compared
     * @return    integer describing the result of the comparison
     */
    public int compareTo(Segment arg) {
        int CompSec = this._Sector < arg._Sector ? -1
                : this._Sector == arg._Sector ? 0 : 1;
        int CompLay = this._Superlayer < arg._Superlayer ? -1
                : this._Superlayer == arg._Superlayer ? 0 : 1;
        int CompId = this._Id < arg._Id ? -1 : this._Id == arg._Id ? 0 : 1;

        int returnVal1 = ((CompId == 0)  ? CompLay : CompId);
        int returnVal2 = ((CompSec == 0) ? returnVal1 : CompSec);

        return returnVal2;
    }

    /**
     * Writes and returns a string containing the segment's information.
     * @return segment's information encoded in a string
     */
    public String getInfo() {
        return "Segment:" +
               "\nID         : " + this.get_Id() +
               "\nSector     : " + this.get_Sector() +
               "\nSuperlayer : " + this.get_Superlayer() +
               "\nSize       : " + this.size();
    }

    /**
     * Writes and returns a string containing all of the segment's data.
     * @return segment's data encoded in a string
     */
    public String getDetailedInfo() {
        return "DC Segment " + this.get_Id() + ":" +
               "\n  Sector             : " + this.get_Sector() +
               "\n  Superlayer         : " + this.get_Superlayer() +
               "\n  Region             : " + this.get_Region() +
               "\n  Residual Sum       : " + this.get_ResiSum() +
               "\n  Time Sum           : " + this.get_TimeSum() +
               "\n  Fit Plane:\n   "       + this.get_fitPlane() +
               "\n  Trajectory         : " + this.get_Trajectory() +
               "\n  Segment End Points : " + this.get_SegmentEndPoints() +
               "\n  Status             : " + this.Status() +
               "\n--------------------------------------------\n";
    }
}
