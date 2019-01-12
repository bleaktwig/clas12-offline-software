package org.jlab.rec.dc.trajectory;

import Jama.*;

/**
 * A StateVec describes a cross measurement in the DC.  It is characterized by a point in the DC
 * tilted coordinate system at each wire plane (i.e. constant z) and by unit tangent vectors in the x and y
 * directions in that coordinate system.  The state vector parameters are the plane index, (x, y, tanTheta_x, tanTheta_y).
 * [The tilted coordinate system is the sector coordinate system rotated by 25 degrees so that the z axis perperdicular to the wire planes]
 * @author ziegler
 *
 */
public class StateVec extends Matrix {

    /**
     * serialVersionUID
     */
    private static final long serialVersionUID = 1874984192960130771L;

    private static int count = 0;
    private int id = -1;

    private int _planeIdx;
    private double z;
    private double b;
    private double h;
    private double _PathLength;

    /**
     * Instantiates a new  vec.
     */
    public StateVec() {
        super(4,1);
        id = count++;
    }

    /**
     * Instantiates a new stateVec
     * @param y the _y
     * @param x the _x
     * @param tanThX the _tanThetaX
     * @param tanThY the _tanThetaY
     */
    public StateVec(double x, double y, double tanThX, double tanThY) {
        super(4,1);
        set(0,0,x);
        set(1,0,y);
        set(2,0,tanThX);
        set(3,0,tanThY);
        id = count++;
    }

    /**
     * Instantiates a new stateVec
     * @param id     the id
     * @param y      the _y
     * @param x      the _x
     * @param tanThX the _tanThetaX
     * @param tanThY the _tanThetaY
     */
    public StateVec(int id, double x, double y, double tanThX, double tanThY) {
        super(4,1);
        set(0,0,x);
        set(1,0,y);
        set(2,0,tanThX);
        set(3,0,tanThY);
        this.id = id;
    }


    /**
     * Instantiates a new stateVec
     *
     * @param v the v
     */
    public StateVec(StateVec v) {
        super(4,1);
        set(0,0,v.x());
        set(1,0,v.y());
        set(2,0,v.tanThetaX());
        set(3,0,v.tanThetaY());
        id = count++;
    }


    /**
     * Instantiates a new  StateVec.
     *
     * @param m the m
     */
    private StateVec(Matrix m) { //needed since Jama.Matrix cannot be casted into StateVec
        super(4,1);
        set(0,0,m.get(0, 0));
        set(1,0,m.get(1, 0));
        set(2,0,m.get(2, 0));
        set(3,0,m.get(3, 0));
        id = count++;
    }

    public int getId() {return id;}
    public void setId(int id) {this.id = id;}

    public void set(StateVec V) {
            set(0,0,V.x());
            set(1,0,V.y());
            set(2,0,V.tanThetaX());
            set(3,0,V.tanThetaY());
    }

    /**
     * Sets the stateVec
     *
     * @param y the _y
     * @param x the _x
     * @param tanThX the _tanThetaX
     * @param tanThY the _tanThetaY
     */
    public void set(double x, double y, double tanThX, double tanThY) {
            set(0,0,x);
            set(1,0,y);
            set(2,0,tanThX);
            set(3,0,tanThY);
    }

    /** Description of x(). */
    public double x() {return(get(0,0));}

    /** Description of y(). */
    public double y() {return(get(1,0));}

    /** Description of tanThetaX(). */
    public double tanThetaX() {return(get(2,0));}

    /** Description of tanThetaY(). */
    public double tanThetaY() {return(get(3,0));}

    public double getPathLength() {return _PathLength;}
    public void setPathLength(double _PathLength) {this._PathLength = _PathLength;}

    /** @return the wire plane index in the series of planes used in the trajectory */
    public int getPlaneIdx() {return _planeIdx;}
    /** Sets the wire plane index in the series of planes used in the trajectory */
    public void setPlaneIdx(int _planeIdx) {this._planeIdx = _planeIdx;}

    public double getZ() {return z;}
    public void setZ(double z) {this.z = z;}

    public double getB() {return b;}
    public void setB(double b) {this.b = b;}

    // KF projector
    public double getProjector() {return h;}
    public void setProjector(double h) {this.h = h;}

    public String getDetailedInfo() {
        return "    StateVec " + this.getId() + ":" +
               "\n      x           : " + this.x() +
               "\n      y           : " + this.y() +
               "\n      tThX        : " + this.tanThetaX() +
               "\n      tThY        : " + this.tanThetaY() +
               "\n      _PathLength : " + this.getPathLength() +
               "\n      _planeIdx   : " + this.getPlaneIdx() +
               "\n      z           : " + this.getZ() +
               "\n      b           : " + this.getB();
    }

    public void printInfo() {
            System.out.println("StateVec [ "+this.x()+", "+this.y()+", "+this.tanThetaX()+", "+this.tanThetaY()+" ] ");
    }

}
