package org.jlab.rec.dc.track;

import java.util.ArrayList;
import java.util.List;
import Jama.Matrix;

import org.jlab.geom.prim.Point3D;
import org.jlab.geom.prim.Vector3D;

import org.jlab.rec.dc.segment.Segment;
import org.jlab.rec.dc.trajectory.StateVec;
import org.jlab.rec.dc.trajectory.Trajectory;

/**
 * A class representing track candidates in the DC. A track has a trajectory represented by an
 * ensemble of geometrical state vectors along its path, a charge and a momentum.
 *
 * @author ziegler
 */
public class Track extends Trajectory implements Comparable<Track>{

    private static final long serialVersionUID = 1763744434903318419L;

    private int _Q;
    private double _P;

    /** Kalman fit covariance matrix. */
    private Matrix _CovMat;

    /** A point along track between the last layer of region 3 and the TOF panel 1b. */
    private Point3D _Region3CrossPoint;

    /** Unit direction vector of the track at a point along the track's trajectory between the last
     *  layer of region 3 and the TOF panel 1b. */
    private Point3D _Region3CrossDir;

    /** Point along the track's trajectory at the HTCC surface sphere. */
    private Point3D _Region1CrossPoint;

    /** Unit direction vector of the track at a point along track trajectory at the HTCC surface
     *  sphere. */
    private Point3D _Region1CrossDir;

    /** Track position at region 1. */
    private Point3D _Region1TrackX;

    /** Track momentum at region 1. */
    private Point3D _Region1TrackP;

    /** The state vector in the tilted sector coordinate system at the mid-plane between the 2
     *  superlayers in region 1. */
    private StateVec _StateVecAtReg1MiddlePlane;
    private int _Id = -1;

    /** Total pathlength of track from vertex to last reference point on trajectory (cm). */
    private double _totPathLen;

    /** Track vertex position at the distance of closest approach to the beam axis (0,0). */
    private Point3D _trakOrig;

    /** Track 3-momentum at the distance of closest approach to the beam axis (0,0). */
    private Vector3D _pOrig;
    private Point3D _Vtx0_TiltedCS;
    private Vector3D _pAtOrig_TiltedCS;

    /** string to indicate if the stage of tracking is hit-based or track-based. */
    private String _trking;

    /** Kalman fit number of degrees of freedom (NDF) */
    private int _FitNDF;

    /** Kalman fit chi^2. */
    private double _fitChisq;
    public boolean fit_Successful;
    private int _missingSuperlayer;

    /** Fit convergence status:
     *  0 : Fit's convergence status OK
     *  1 : The fit exits before converging.
     */
    private int _fitConvergenceStatus;

    private int _Status = 0;

    /** Associated track for Hit-Based tracking. */
    private Track _AssociatedHBTrack;

    /** List of HB Segments, used for Time-Based tracking. */
    private List<Segment> _ListOfHBSegments = new ArrayList<Segment>();

    /** empty constructor. */
    public Track() {}

    public int get_MissingSuperlayer() {return _missingSuperlayer;}
    public void set_MissingSuperlayer(int missingSuperlayer) {
        this._missingSuperlayer = missingSuperlayer;
    }

    public int get_Status() {return _Status;}
    public void set_Status(int _Status) {this._Status = _Status;}

    public int get_Id() {return _Id;}
    public void set_Id(int _Id) {this._Id = _Id;}

    public int get_Q() {return _Q;}
    public void set_Q(int _Q) {this._Q = _Q;}

    public double get_P() {return _P;}
    public void set_P(double _P) {this._P = _P;}

    public Point3D get_PostRegion3CrossPoint() {return _Region3CrossPoint;}
    public void set_PostRegion3CrossPoint(Point3D point) {_Region3CrossPoint = point;}

    public Point3D get_PostRegion3CrossDir() {return _Region3CrossDir;}
    public void set_PostRegion3CrossDir(Point3D dir) {_Region3CrossDir = dir;}

    public Point3D get_PreRegion1CrossPoint() {return _Region1CrossPoint;}
    public void set_PreRegion1CrossPoint(Point3D point) {_Region1CrossPoint = point;}

    public Point3D get_PreRegion1CrossDir() {return _Region1CrossDir;}
    public void set_PreRegion1CrossDir(Point3D dir) {_Region1CrossDir = dir;}

    public Point3D get_Region1TrackX() {return _Region1TrackX;}
    public void set_Region1TrackX(Point3D _Region1TrackX) {this._Region1TrackX = _Region1TrackX;}

    public Point3D get_Region1TrackP() {return _Region1TrackP;}
    public void set_Region1TrackP(Point3D _Region1TrackP) {this._Region1TrackP = _Region1TrackP;}

    public double get_TotPathLen() {return _totPathLen;}
    public void set_TotPathLen(double totPathLen) {_totPathLen = totPathLen;}

    public Point3D get_Vtx0() {return _trakOrig;}
    public void set_Vtx0(Point3D trakOrig) {_trakOrig = trakOrig;}

    public Vector3D get_pAtOrig() {return _pOrig;}
    public void set_pAtOrig(Vector3D pOrig) {_pOrig = pOrig;}

    public String get_TrackingInfoString() {return _trking;}
    public void set_TrackingInfoString(String trking) {_trking = trking;}

    public int get_FitNDF() {return _FitNDF;}
    public void set_FitNDF(int _FitNDF) {this._FitNDF = _FitNDF;}

    public double get_FitChi2() {return _fitChisq;}
    public void set_FitChi2(double fitChisq) {_fitChisq = fitChisq;}

    public Matrix get_CovMat() {return _CovMat;}
    public void set_CovMat(Matrix _CovMat) {this._CovMat = _CovMat;}

    public int get_FitConvergenceStatus() {return _fitConvergenceStatus;}
    public void set_FitConvergenceStatus(int fitConvergenceStatus) {
        this._fitConvergenceStatus = fitConvergenceStatus;
    }

    public StateVec get_StateVecAtReg1MiddlePlane() {return _StateVecAtReg1MiddlePlane;}
    public void set_StateVecAtReg1MiddlePlane(StateVec _StateVecAtReg1MiddlePlane) {
        this._StateVecAtReg1MiddlePlane = _StateVecAtReg1MiddlePlane;
    }

    public Track get_AssociatedHBTrack() {return _AssociatedHBTrack;}
    public void set_AssociatedHBTrack(Track _trk) {_AssociatedHBTrack = _trk;}

    public List<Segment> get_ListOfHBSegments() {return _ListOfHBSegments;}
    public void set_ListOfHBSegments(List<Segment> _listOfHBSegments) {
        this._ListOfHBSegments = _listOfHBSegments;
    }


    /** Prints basic track info */
    public void printInfo() {
        System.out.println("Track " + this._Id + " Q= " + this._Q + " P= " + this._P);
    }

    /** Prints every variable in the track */
    public void printDetailedInfo() {
        System.out.println("Track " + this._Id + ":");
        System.out.println("  _MissingSuperlayer         : " + get_MissingSuperlayer());
        System.out.println("  _Status                    : " + get_Status());
        System.out.println("  _Q                         : " + get_Q());
        System.out.println("  _P                         : " + get_P());
        System.out.println("  _Region3CrossPoint         : " + get_PostRegion3CrossPoint()); // POINT/DIR
        System.out.println("  _Region3CrossDir           : " + get_PostRegion3CrossDir()); // POINT/DIR
        System.out.println("  _Region1CrossPoint         : " + get_PreRegion1CrossPoint()); // POINT/DIR
        System.out.println("  _Region1CrossDir           : " + get_PreRegion1CrossDir()); // POINT/DIR
        System.out.println("  _Region1TrackX             : " + get_Region1TrackX()); // POINT/DIR
        System.out.println("  _Region1TrackP             : " + get_Region1TrackP()); // POINT/DIR
        System.out.println("  _TotPathLen                : " + get_TotPathLen());
        System.out.println("  _Vtx0                      : " + get_Vtx0()); // POINT/DIR
        System.out.println("  _pOrig                     : " + get_pAtOrig()); // VECTOR
        System.out.println("  _TrackingInfoString        : " + get_TrackingInfoString());
        System.out.println("  _FitNDF                    : " + get_FitNDF());
        System.out.println("  _fitChisq                  : " + get_FitChi2());
        System.out.println("  _CovMat                    : "); // MATRIX
        String CovMatStr = "";
        for (int ii = 0; ii < get_CovMat().getRowDimension(); ++ii) {
            CovMatStr += "    ";
            for (int jj = 0; jj < get_CovMat().getColumnDimension(); ++jj) {
                CovMatStr += get_CovMat().get(ii, jj);
                CovMatStr += " ";
            }
            CovMatStr += "\n";
        }
        System.out.println(CovMatStr);
        System.out.println("  _FitConvergenceStatus      : " + get_FitConvergenceStatus());
        System.out.println("  _StateVecAtReg1MiddlePlane : " + get_StateVecAtReg1MiddlePlane()); // STATEVEC
        if (get_StateVecAtReg1MiddlePlane() != null) {
            System.out.println(get_StateVecAtReg1MiddlePlane().getDetailedInfo());
        }
        else {
            System.out.println("  _StateVecAtReg1MiddlePlane : null");
        }
        if (get_AssociatedHBTrack() != null)
            System.out.println("  _AssociatedHBTrack         : " + get_AssociatedHBTrack()._Id);
        else
            System.out.println("  _AssociatedHBTrack         : null");
        System.out.println("  list of segments:");
        if (get_ListOfHBSegments().size() == 0) System.out.println("    empty.");
        else {
            for (Segment segment : get_ListOfHBSegments())
                System.out.println("    Segment " + segment.get_Id());
        }
    }

    // NOTE: Lacks explanation of the comparison criteria.
    /**
     * Compares the track to another one.
     * @param  arg the track to which this one is compared
     * @return     integer given by comparison criteria
     */
    @Override
    public int compareTo(Track arg) {
        int idtrkSeg1 = this.get(0).get_Segment1().get(0).get_AssociatedHBTrackID();
        int idtrkSeg1a = arg.get(0).get_Segment1().get(0).get_AssociatedHBTrackID();
        int idtrkSeg2 = this.get(0).get_Segment2().get(0).get_AssociatedHBTrackID();
        int idtrkSeg2a = arg.get(0).get_Segment2().get(0).get_AssociatedHBTrackID();

        int idtrkSeg3 = this.get(1).get_Segment1().get(0).get_AssociatedHBTrackID();
        int idtrkSeg3a = arg.get(1).get_Segment1().get(0).get_AssociatedHBTrackID();
        int idtrkSeg4 = this.get(1).get_Segment2().get(0).get_AssociatedHBTrackID();
        int idtrkSeg4a = arg.get(1).get_Segment2().get(0).get_AssociatedHBTrackID();

        int idtrkSeg5 = this.get(2).get_Segment1().get(0).get_AssociatedHBTrackID();
        int idtrkSeg5a = arg.get(2).get_Segment1().get(0).get_AssociatedHBTrackID();
        int idtrkSeg6 = this.get(2).get_Segment2().get(0).get_AssociatedHBTrackID();
        int idtrkSeg6a = arg.get(2).get_Segment2().get(0).get_AssociatedHBTrackID();


        int return_val1 = idtrkSeg1 < idtrkSeg1a ? -1 : idtrkSeg1 == idtrkSeg1a ? 0 : 1;
        int return_val2 = idtrkSeg2 < idtrkSeg2a ? -1 : idtrkSeg2 == idtrkSeg2a ? 0 : 1;
        int return_val3 = idtrkSeg3 < idtrkSeg3a ? -1 : idtrkSeg3 == idtrkSeg3a ? 0 : 1;
        int return_val4 = idtrkSeg4 < idtrkSeg4a ? -1 : idtrkSeg4 == idtrkSeg4a ? 0 : 1;
        int return_val5 = idtrkSeg5 < idtrkSeg5a ? -1 : idtrkSeg5 == idtrkSeg5a ? 0 : 1;
        int return_val6 = idtrkSeg6 < idtrkSeg6a ? -1 : idtrkSeg6 == idtrkSeg6a ? 0 : 1;

        int return_val_a1 = ((return_val1 == 0) ? return_val2   : return_val1);
        int return_val_a2 = ((return_val2 == 0) ? return_val_a1 : return_val2);
        int return_val_a3 = ((return_val3 == 0) ? return_val_a2 : return_val3);
        int return_val_a4 = ((return_val4 == 0) ? return_val_a3 : return_val4);
        int return_val_a5 = ((return_val5 == 0) ? return_val_a4 : return_val5);
        int return_val_a6 = ((return_val6 == 0) ? return_val_a5 : return_val6);

        int returnSec = this.get_Sector() < arg.get_Sector() ? -1
                      : this.get_Sector() == arg.get_Sector() ? 0 : 1;

        return ((returnSec == 0) ? return_val_a6 : returnSec);
    }
}
