package org.jlab.rec.dc.hit;

import eu.mihosoft.vrl.v3d.Vector3d;

import org.jlab.geom.prim.Point3D;
import org.jlab.utils.groups.IndexedTable;
import org.jlab.detector.geant4.v2.DCGeant4Factory;

import org.jlab.rec.dc.Constants;
import org.jlab.rec.dc.timetodistance.TimeToDistanceEstimator;
import org.jlab.rec.dc.trajectory.StateVec;

/**
 * A hit that was used in a fitted cluster. It extends the Hit class and contains local and sector
 * coordinate information at the MidPlane. An estimate for the Left-right Ambiguity is assigned
 * based on the linear fit to the wire position residual.
 *
 * @author ziegler
 *
 */
public class FittedHit extends Hit implements Comparable<Hit> {

    private double _X;                // The hit's x-position at the mid-plane (y = 0) in the tilted
                                      //     sector's coordinate system.
    private double _XMP;              // X at the MidPlane in sector coordinates system.
    private double _Z;                // The hit's z-position at the mid-plane (y = 0) in the tilted
                                      //     sector's coordinate system.
    private double _lX;			      // X in local coordinate system used in hit-based fit to
                                      //     cluster line, used in the cluster-finding algorithm to
                                      //     fit the hit-based wire positions, from 1 to 6.
    private double _lY;			      // Y in local coordinate system used in hit-based fit to
                                      //     cluster line, used in the cluster-finding algorithm to
                                      //     fit the hit-based wire positions.
    private double _Residual;		  // Cluster line to the wire position's residual, or the
                                      // residual from the fit to the wire positions in the
                                      // superlayer.
    private double _TimeResidual = 0; // Cluster line to the wire position's time-residual, or
                                      //     |fit| - |y| from the fit to the wire positions in the
                                      //     superlayer.
    private int _LeftRightAmb;		  // Left-Right Ambiguity value. An integer representative of
                                      //     estimate of the left-right ambiguity obtained from
                                      //     pattern recognition:
                                      //     -1: the track went to the left of the wire.
                                      //      0: the left-right ambiguity could not be resolved.
                                      //     +1: the track went to the right of the wire.
    private double _QualityFac;       // A quality factor representative of the quality of the fit
                                      //     to the hit.
    private int _TrkgStatus = -1;	  // An integer representative of the pattern recognition and
                                      //     subsequent KF fit for the hit.
                                      //     -1: the hit has not yet been fit and is the input of
                                      //         hit-based tracking, having a well-defined time-to-
                                      //         distance value.
                                      //      0: the hit has been successfully involved in hit-based
                                      //         tracking and has a wel-defined time-to-distance.
                                      //     +1: the hit has been succesfully involved in track
                                      //         fitting.

    private double _ClusFitDoca = -1;   // Doca to cluster fit line in cm.
    private double _TrkFitDoca = -1;    // Doca to track trajectory at hit layer plane in cm.
    private double _TimeToDistance = 0; // the calculated distance in cm from the time in ns.
    private double _Beta = 1.0;         // NOTE: Needs description

    private StateVec _AssociatedStateVec; // State vector (x, y, tx, ty, q/p) associated with the
                                          //     hit.

    private double _Doca;	 // Reconstructed doca, for now it is using the linear parametrization
                             //     that is in gemc. Measured in cm.
    // private double _DocaErr; // Error on doca.
    private double _B;		 // B-field intensity at hit location along wire, measured in T.
    private int _Id;         // The ID corresponds to the hit index in the EvIO column.
    public int _lr;          // NOTE: Needs description

    public boolean RemoveFlag = false;
    private int _AssociatedClusterID = -1; // ID of the cluster associated to the fitted hit.
    private int _AssociatedHBTrackID = -1; // ID of the HB track associated to the fitted hit.
    private int _AssociatedTBTrackID = -1; // ID of the TB track associated to the fitted hit.

    // intersection of cross direction line with the hit wire (TCS)
    private Point3D CrossDirIntersWire;        // NOTE: Needs description
    private double _SignalPropagAlongWire;     // NOTE: Needs description
    private double _SignalPropagTimeAlongWire; // NOTE: Needs description
    private double _SignalTimeOfFlight;        // Signal time of flight to the track doca of the hit
                                               // in ns.

    private double _T0;      // T0 calibration constant in ns
    private double _tFlight; // Flight time to the track's closest point to the hit wire in ns
    private double _tProp;   // Propagation time along the wire in ns
    private double _tStart;  // The event start time in ns from EB bank.
    private double _Time;    // Time = TDC - tFlight - tProp - T0 - TStart in ns

    private boolean _OutOfTimeFlag; // Boolean flag for identifying out of time hits.
    private double _WireLength;     // Length of the wire.
    private double _WireMaxSag;     // NOTE: Needs description
    private double _TrkResid = 999; // NOTE: Needs description

    private double _deltatime_beta; // NOTE: Needs description

    public FittedHit(int sector, int superlayer, int layer, int wire, int TDC, int id) {
        super(sector, superlayer, layer, wire, TDC, id);

        this.set_lX(layer);
        this.set_lY(layer, wire);
    }

    public double getB() {return _B;}
    public void setB(double _B) {this._B = _B;}

    public int get_Id() {return _Id;}
    public void set_Id(int _Id) {this._Id = _Id;}

    public double get_Doca() {return _Doca;}
    public void set_Doca(double _Doca) {this._Doca = _Doca;}

    public double get_lX() {return _lX;}
    public void set_lX(double layerValue) {this._lX = layerValue;}

    public double get_lY() {return _lY;}

    /** @param lY explicit definition for lY */
    public void set_lY(double _lY) {this._lY = _lY;}

    /**
     * Calculates the center of the cell utilizing the layer and the wire in the local superlayer
     * coordinate system and sets it.
     * @param layer layer number from 1 to 6
     * @param wire  wire number from 1 to 112 sets the center of the cell as a function of wire
     *              number in the local superlayer coordinate system.
     */
    public void set_lY(int layer, int wire) {this._lY = this.calcLocY(layer, wire);}

    public double get_TimeResidual() {return _TimeResidual;}
    public void set_TimeResidual(double _TimeResidual) {this._TimeResidual = _TimeResidual;}

    public double get_Residual() {return _Residual;}
    public void set_Residual(double _Residual) {this._Residual = _Residual;}

    public int get_LeftRightAmb() {return _LeftRightAmb;}
    public void set_LeftRightAmb(int leftRightAmb) {this._LeftRightAmb = leftRightAmb;}

    public double get_QualityFac() {return _QualityFac;}
    public void set_QualityFac(double _QualityFac) {this._QualityFac = _QualityFac;}

    public int get_TrkgStatus() {return _TrkgStatus;}
    public void set_TrkgStatus(int trkgStatus) {_TrkgStatus = trkgStatus;}

    public double get_TimeToDistance() {return _TimeToDistance;}

    /**
     * Sets the calculated distance (in cm) from the time (in ns).
     * @param cosTrkAngle NOTE: Missing description
     * @param B           NOTE: Missing description
     * @param tab         NOTE: Missing description
     * @param tde         NOTE: Missing description
     */
    public void set_TimeToDistance(double cosTrkAngle,
                                   double B,
                                   IndexedTable tab,
                                   TimeToDistanceEstimator tde) {

        double distance = 0;
        int slIdx  = this.get_Superlayer() - 1;
        int secIdx = this.get_Sector() - 1;

        if (_TrkgStatus != -1 && this.get_Time() > 0) {
            double alpha  = Math.acos(cosTrkAngle);
            double ralpha = this.reducedAngle(alpha);
            double beta   = this.get_Beta();
            double x      = this.get_ClusFitDoca();
            double deltatime_beta = 0;

            if (x != -1) {
                double V_0 = Constants.V0AVERAGED;
                deltatime_beta = (Math.sqrt(x * x +
                        (tab.getDoubleValue("distbeta", this.get_Sector(),
                                            this.get_Superlayer(),0) * beta * beta) *
                        (tab.getDoubleValue("distbeta", this.get_Sector(),
                                            this.get_Superlayer(),0) * beta * beta)) - x
                                 ) / V_0;
            }
            this.set_DeltaTimeBeta(deltatime_beta);
            double correctedTime = (this.get_Time() - deltatime_beta);
            if (correctedTime <= 0) correctedTime = 0.01;

            distance = tde.interpolateOnGrid(B, Math.toDegrees(ralpha), correctedTime, secIdx, slIdx);
        }

        this.set_Doca(distance);
        this._TimeToDistance = distance;
    }

    /**
     * NOTE: Missing description
     * @param cellSize the cell size in cm
     */
    public void fix_TimeToDistance(double cellSize) {this._TimeToDistance = cellSize;}

    public StateVec getAssociatedStateVec() {return _AssociatedStateVec;}
    public void setAssociatedStateVec(StateVec _AssociatedStateVec) {
        this._AssociatedStateVec = _AssociatedStateVec;
    }

    public double get_ClusFitDoca() {return _ClusFitDoca;}
    public void set_ClusFitDoca(double _ClusFitDoca) {this._ClusFitDoca = _ClusFitDoca;}

    public double get_TrkFitDoca() {return _TrkFitDoca;}
    public void set_TrkFitDoca(double _TrkFitDoca) {this._TrkFitDoca = _TrkFitDoca;}

    public double get_X() {return _X;}
    public void set_X(double _X) {this._X = _X;}

    public double get_XMP() {return _XMP;}
    public void set_XMP(double _XMP) {this._XMP = _XMP;}

    public double get_Z() {return _Z;}
    public void set_Z(double _Z) {this._Z = _Z;}

    public double get_WireLength() {return _WireLength;}
    public void set_WireLength(double _WireLength) {this._WireLength = _WireLength;}

    public double get_WireMaxSag() {return _WireMaxSag;}
    public void set_WireMaxSag(double _WireMaxSag) {this._WireMaxSag = _WireMaxSag;}

    public double get_TrkResid() {return _TrkResid;}
    public void set_TrkResid(double _TrkResid) {this._TrkResid = _TrkResid;}

    public int get_AssociatedClusterID() {return _AssociatedClusterID;}
    public void set_AssociatedClusterID(int _AssociatedClusterID) {
        this._AssociatedClusterID = _AssociatedClusterID;
    }

    public int get_AssociatedHBTrackID() {return _AssociatedHBTrackID;}
    public void set_AssociatedHBTrackID(int _id) {_AssociatedHBTrackID = _id;}

    public int get_AssociatedTBTrackID() {return _AssociatedTBTrackID;}
    public void set_AssociatedTBTrackID(int _id) {_AssociatedTBTrackID = _id;}

    public Point3D getCrossDirIntersWire() {return CrossDirIntersWire;}
    public void setCrossDirIntersWire(Point3D CrossDirIntersWire) {
        this.CrossDirIntersWire = CrossDirIntersWire;
    }

    public double get_Beta() {return _Beta;}
    public void set_Beta(double beta) {_Beta = beta;}

    public double getSignalPropagAlongWire() {return _SignalPropagAlongWire;}
    public void setSignalPropagAlongWire(DCGeant4Factory DcDetector) {
        this._SignalPropagAlongWire = this.calc_SignalPropagAlongWire(DcDetector);
    }

    /**
     * Calculates the signal propagation time along the wire in ns.
     * @param DcDetector DC detector geometry
     * @return           signal propagation time along the wire
     */
    public double calc_SignalPropagAlongWire(DCGeant4Factory DcDetector) {

        Vector3d WireEnd;
        if (Constants.STBLOC[this.get_Sector() - 1][this.get_Superlayer() - 1] > 0) {
            WireEnd = DcDetector.getWireRightend(this.get_Sector()-1,  this.get_Superlayer() - 1,
                                                 this.get_Layer() - 1, this.get_Wire() - 1);
        } else {
            WireEnd = DcDetector.getWireLeftend(this.get_Sector()-1,  this.get_Superlayer() - 1,
                                                this.get_Layer() - 1, this.get_Wire() - 1);
        }

        double X = this.getCrossDirIntersWire().x();
        double Y = this.getCrossDirIntersWire().y();

        return Math.sqrt((X - WireEnd.x) * (X - WireEnd.x) + (Y - WireEnd.y) * (Y - WireEnd.y));
    }

    /**
     * Calculates the signal propagation time along the wire in ns given an explicit X and Y.
     * @param X          NOTE: Missing description
     * @param Y          NOTE: Missing description
     * @param DcDetector DC detector geometry
     * @return           signal propagation time along the wire
     */
    public double calc_SignalPropagAlongWire(double X, double Y, DCGeant4Factory DcDetector) {

        Vector3d WireEnd;
        if (Constants.STBLOC[this.get_Sector()-1][this.get_Superlayer()-1] > 0) {
            WireEnd = DcDetector.getWireRightend(this.get_Sector() - 1, this.get_Superlayer() - 1,
                                                 this.get_Layer() - 1,  this.get_Wire() - 1);
        } else {
            WireEnd = DcDetector.getWireLeftend(this.get_Sector() - 1, this.get_Superlayer() - 1,
                                                this.get_Layer() - 1,  this.get_Wire() - 1);
        }

        return Math.sqrt((X - WireEnd.x) * (X - WireEnd.x) + (Y - WireEnd.y) * (Y - WireEnd.y));
    }

    public double getSignalPropagTimeAlongWire() {return _SignalPropagTimeAlongWire;}
    public void setSignalPropagTimeAlongWire(DCGeant4Factory DcDetector) {
        this.setSignalPropagAlongWire(DcDetector);
        this._SignalPropagTimeAlongWire = this._SignalPropagAlongWire / (Constants.SPEEDLIGHT*0.7);
        this._tProp = this._SignalPropagTimeAlongWire;
    }

    public double getSignalTimeOfFlight() {return _SignalTimeOfFlight;}

    /** Sets signal time of flight to the track doca to the hit wire in ns. */
    public void setSignalTimeOfFlight() {
        if (this.get_Beta() > 0 && this.getAssociatedStateVec() != null) {
            this._SignalTimeOfFlight = (this.getAssociatedStateVec().getPathLength()) /
                                       (Constants.SPEEDLIGHT * this.get_Beta());
        }
        this._tFlight = this._SignalTimeOfFlight;
    }

    public double getTStart() {return _tStart;}
    public void setTStart(double tStart) {this._tStart = tStart;}

    public double getT0() {return _T0;}
    public void setT0(double T0) {this._T0 = T0;}

    public double getTFlight() {return _tFlight;}
    public void setTFlight(double tFlight) {this._tFlight = tFlight;}

    public double getTProp() {return _tProp;}
    public void setTProp(double tProp) {this._tProp = tProp;}

    public double get_Time() {return _Time;}
    public void set_Time(double _Time) {this._Time = _Time;}

    public boolean get_OutOfTimeFlag() {return _OutOfTimeFlag;}
    public void set_OutOfTimeFlag(boolean b) {_OutOfTimeFlag = b;}

    public double get_DeltaTimeBeta() {return _deltatime_beta ;}
    public void set_DeltaTimeBeta(double deltatime_beta) {_deltatime_beta = deltatime_beta;}

    /**
     * Calculates the approximate uncertainty on the hit's positions using the inverse of the gemc
     * smearing function.
     * @param B          NOTE: Missing description
     * @param constants0 NOTE: Missing description
     * @param constants1 NOTE: Missing description
     * @param tde        NOTE: Missing description
     * @return           the approximate uncertainty on the hit position
     */
    public double get_PosErr(double B,
                             IndexedTable constants0,
                             IndexedTable constants1,
                             TimeToDistanceEstimator tde) {

        double err = this.get_DocaErr();

        if (this._TrkgStatus != -1) {
            if (this.get_TimeToDistance() == 0) {
                // If the time-to-dist is not set, set it
                set_TimeToDistance(1.0, B, constants1, tde);
            }

            err = Constants.CELLRESOL; // Default
            double x = this.get_Doca() / this.get_CellSize();
            double p1 = constants0.getDoubleValue("parameter1",
                                                  this.get_Sector(),this.get_Superlayer(), 0);
            double p2 = constants0.getDoubleValue("parameter2",
                                                  this.get_Sector(),this.get_Superlayer(), 0);
            double p3 = constants0.getDoubleValue("parameter3",
                                                  this.get_Sector(),this.get_Superlayer(), 0);
            double p4 = constants0.getDoubleValue("parameter4",
                                                  this.get_Sector(),this.get_Superlayer(), 0);
            double scale = constants0.getDoubleValue("scale",
                                                     this.get_Sector(),this.get_Superlayer(), 0);

            // Gives a reasonable approximation to the measured CLAS resolution
            //     in cm! --> scale by 0.1
            err = (p1 + p2 / ((p3 + x) * (p3 + x)) + p4 * Math.pow(x, 8)) * scale * 0.1;
        }

        return err;
    }

    /**
     * Reduces a given angle between the range of 0 and 30 degrees.
     * @param alpha the local angle of the track
     * @return      the reduced angle
     */
    double reducedAngle(double alpha) {
        // Math.PI / 3. = 1.0471975511965976
        // Math.PI / 6. = 0.5235987755982988

        double ralpha = Math.abs(alpha);

        while (ralpha > 1.0471975511965976) {
            ralpha -= 1.0471975511965976;
        }

        if (ralpha > 0.5235987755982988) {
            ralpha = 1.0471975511965976 - ralpha;
        }

        return ralpha;
    }

    /**
     * A method to update the hit position information after the fit to the local coordinate
     * system's wire positions.
     * @param DcDetector DC detector geometry
     */
    public void updateHitPosition(DCGeant4Factory DcDetector) {
        double z = DcDetector.getWireMidpoint(this.get_Sector() - 1,
                                              this.get_Superlayer() - 1,
                                              this.get_Layer() - 1,
                                              this.get_Wire() - 1).z;
        double x = this.calc_GeomCorr(DcDetector, 0);

        this.set_X(x);
        this.set_Z(z);
    }

    /**
     * A method to update the hit position information after the fit to the wire positions employing
     * hit-based tracking algorithms has been performed.
     * @param cosTrkAngle NOTE: Missing description
     * @param B           NOTE: Missing description
     * @param tab         NOTE: Missing description
     * @param DcDetector  DC detector geometry
     * @param tde         NOTE: Missing description
     */
    public void updateHitPositionWithTime(double cosTrkAngle,
                                          double B,
                                          IndexedTable tab,
                                          DCGeant4Factory DcDetector,
                                          TimeToDistanceEstimator tde) {

        // System.out.println("[DCKF] TIME: " + this.get_Time());
        if (this.get_Time() > 0) this.set_TimeToDistance(cosTrkAngle, B, tab, tde);

        double z = DcDetector.getWireMidpoint(this.get_Sector() - 1,
                                              this.get_Superlayer() - 1,
                                              this.get_Layer() - 1,
                                              this.get_Wire() - 1).z;
        double x = this.calc_GeomCorr(DcDetector, 0);
        double MPCorr = 1;
        if (cosTrkAngle > 0.8 && cosTrkAngle <= 1) MPCorr = cosTrkAngle;

        // Math.cos(Math.toRadians(6.)) = 0.9945218953682733
        this.set_X(x + this.get_LeftRightAmb() *
                       (this.get_TimeToDistance() / MPCorr) / 0.9945218953682733);
        this.set_Z(z);
    }

    public double XatY(DCGeant4Factory DcDetector, double y) {
        return this.calc_GeomCorr(DcDetector, y) +
               this.get_LeftRightAmb() * (this.get_TimeToDistance());
    }

    private double calc_GeomCorr(DCGeant4Factory DcDetector, double y) {

        double xL = DcDetector.getWireLeftend(this.get_Sector() - 1, this.get_Superlayer() - 1,
                                              this.get_Layer() - 1,  this.get_Wire() - 1).x;
        double xR = DcDetector.getWireRightend(this.get_Sector() - 1, this.get_Superlayer() - 1,
                                               this.get_Layer() - 1,  this.get_Wire() - 1).x;
        double yL = DcDetector.getWireLeftend(this.get_Sector() - 1, this.get_Superlayer() - 1,
                                              this.get_Layer() - 1,  this.get_Wire() - 1).y;
        double yR = DcDetector.getWireRightend(this.get_Sector() - 1, this.get_Superlayer() - 1,
                                               this.get_Layer() - 1, this.get_Wire() - 1).y;

        double DL = Constants.MAXENDPLTDEFLEC[this.get_Region() - 1][this.get_Sector() - 1][0];
        double DR = Constants.MAXENDPLTDEFLEC[this.get_Region() - 1][this.get_Sector() - 1][1];

        double wire = this.get_Wire();
        double normW = (double) wire/112.;

        xL -= Constants.getWIREDIST()*DL*(normW - 3*normW*normW*normW + 2*normW*normW*normW*normW);
        xR -= Constants.getWIREDIST()*DR*(normW - 3*normW*normW*normW + 2*normW*normW*normW*normW);

        double x = xR - (yR - y)*((xR - xL)/(yR - yL));
        double wireLen = Math.sqrt((xL - xR)*(xL - xR) + (yL - yR)*(yL - yR));

        int sector = this.get_Sector();
        int A = 0;
        switch (sector) {
            case (1): A = 0;
                      break;
            case (2): A = -1;
                      break;
            case (3): A = -1;
                      break;
            case (4): A = 0;
                      break;
            case (5): A = 1;
                      break;
            case (6): A = 1;
                      break;
            default:  throw new RuntimeException("invalid sector");
        }

        int region = this.get_Region();
        double C = 0;
        double ConvFac = 1000000;
        switch (region) {
            case (1): C = 2.0/ConvFac;
                      break;
            case (2): C = 4.95/ConvFac;
                      break;
            case (3): if      (wire < 69) C = 12.5/ConvFac;
                      else if (wire < 92) C = 7.49/ConvFac;
                      else                C = 5.98/ConvFac;
                      break;
            default:  throw new RuntimeException("invalid region");
        }

        // Math.cos(Math.toRadians(25.)) = 0.9063077870366499
        // Math.cos(Math.toRadians(30.)) = 0.8660254037844387
        double MaxSag = Constants.getWIREDIST() * A * C * wire * wire *
                        0.9063077870366499 * 0.8660254037844387;

        double delta_x = MaxSag * (1. - y/(0.5*wireLen)) * (1. - y/(0.5*wireLen));

        x += delta_x;

        this.set_WireLength(wireLen);
        this.set_WireMaxSag(MaxSag);

        return x;
    }

    /**
     * Compares to hits based on basic descriptors, returning true if the hits are the same and
     * false otherwise.
     * @param otherHit the other hit
     * @return         a boolean describing the result of the comparison
     */
    public boolean isSameAs(FittedHit otherHit) {
        boolean cmp = false;
        if (this.get_Time() == otherHit.get_Time()
                && this.get_Sector() == otherHit.get_Sector()
                && this.get_Superlayer() == otherHit.get_Superlayer()
                && this.get_Layer() == otherHit.get_Layer()
                && this.get_Wire() == otherHit.get_Wire()) {
            cmp = true;
        }
        return cmp;
    }

    /**
     * Compares hits to sort them by layer number.
     * @param otherHit the other hit
     * @return         an int used to sort a collection of hits
     */
    public int compareTo(FittedHit otherHit) {
        if (this.get_Layer() > otherHit.get_Layer()) return 1;
        else                                         return -1;
    }

    /**
     * Writes and returns a string containing the fitted hit's information.
     * @return hit's information encoded in a string
     */
    public String getInfo() {
        return "DC Fitted Hit:" +
               "\nID          : " + this.get_Id() +
               "\nSector      : " + this.get_Sector() +
               "\nSuperlayer  : " + this.get_Superlayer() +
               "\nLayer       : " + this.get_Layer() +
               "\nWire        : " + this.get_Wire() +
               "\nTDC         : " + this.get_TDC() +
               "\nTime        : " + this.get_Time() +
               "\nLR          : " + this.get_LeftRightAmb() +
               "\ndoca        : " + this.get_TimeToDistance() +
               "\n+/-         : " + this.get_DocaErr() +
               "\nupdated pos : " + this._X +
               "\nclus        : " + this._AssociatedClusterID;
    }

    /**
     * Writes and returns a string containing all of the fitted hit's data.
     * @return fitted hit's data encoded in a string
     */
    public String getDetailedInfo() {
        String returnStr = "DC Fitted Hit " + this.get_Id() + ":" +
                           "\n  Hit Data:" +
                           "\n    Sector                        : " + this.get_Sector() +
                           "\n    Superlayer                    : " + this.get_Superlayer() +
                           "\n    Layer                         : " + this.get_Layer() +
                           "\n    Wire                          : " + this.get_Wire() +
                           "\n    TDC                           : " + this.get_TDC() +
                           "\n    Cell size                     : " + this.get_CellSize() +
                           "\n    Doca err                      : " + this.get_DocaErr() +
                           "\n  Fitted Hit Data:" +
                           "\n    B                             : " + this.getB() +
                           "\n    Doca (_Doca)                  : " + this.get_Doca() +
                           "\n    Doca (_TimeToDistance)        : " + this.get_TimeToDistance() +
                           "\n    Doca Error                    : " + this.get_DocaErr() +
                           "\n    _lX                           : " + this.get_lX() +
                           "\n    _lY                           : " + this.get_lY() +
                           "\n    Time Residual                 : " + this.get_TimeResidual() +
                           "\n    Residual                      : " + this.get_Residual() +
                           "\n    Left Right Ambiguity          : " + this.get_LeftRightAmb() +
                           "\n    Quality Fac                   : " + this.get_QualityFac() +
                           "\n    Tracking Status               : " + this.get_TrkgStatus() +
                           "\n    Time to Distance              : " + this.get_TimeToDistance();
        if (this.getAssociatedStateVec() != null)
            returnStr += "\n    Associated State Vector         :\n" +
                         this.getAssociatedStateVec().getDetailedInfo();
        else
            returnStr += "\n    Associated State Vector         : null";
        returnStr += "\n    Cluster Fit Doca              : " + this.get_ClusFitDoca() +
                     "\n    Tracking Fit Doca             : " + this.get_TrkFitDoca() +
                     "\n    Updated Position (_X)         : " + this.get_X() +
                     "\n    _XMP                          : " + this.get_XMP() +
                     "\n    _Z                            : " + this.get_Z() +
                     "\n    Wire Length                   : " + this.get_WireLength() +
                     "\n    _WireMaxSag                   : " + this.get_WireMaxSag() +
                     "\n    Tracking Residual             : " + this.get_TrkResid() +
                     "\n    Associated Cluster Id         : " + this.get_AssociatedClusterID() +
                     "\n    Associated HB Track Id        : " + this.get_AssociatedHBTrackID() +
                     "\n    Associated TB Track Id        : " + this.get_AssociatedTBTrackID() +
                     "\n    Cross Dir Inters Wire         : " + this.getCrossDirIntersWire() +
                     "\n    Beta                          : " + this.get_Beta() +
                     "\n    Signal Propag Along Wire      : " + this.getSignalPropagAlongWire() +
                     "\n    Signal Propag Time Along Wire : " + this.getSignalPropagTimeAlongWire() +
                     "\n    Signal Time of Flight         : " + this.getSignalTimeOfFlight() +
                     "\n    T Start                       : " + this.getTStart() +
                     "\n    T0                            : " + this.getT0() +
                     "\n    TFlight                       : " + this.getTFlight() +
                     "\n    TProp                         : " + this.getTProp() +
                     "\n    Time (_Time)                  : " + this.get_Time() +
                     "\n    Out of Time Flag              : " + this.get_OutOfTimeFlag() +
                     "\n    Delta Time Beta               : " + this.get_DeltaTimeBeta() +
                     // "\n    Position Error                : " + this.get_PosErr() +
                     "\n----------------------\n";
        return returnStr;
    }
}
