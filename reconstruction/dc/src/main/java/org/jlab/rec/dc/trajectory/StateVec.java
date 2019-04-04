package org.jlab.rec.dc.trajectory;

import Jama.Matrix;

/**
 * A StateVec describes a cross measurement in the DC. It is characterized by a point in the DC
 * tilted coordinate system at each wire plane (i.e. constant z) and by unit tangent vectors in the
 * x and y directions in that coordinate system. The state vector parameters are the plane index
 * (x, y, tanTheta_x, tanTheta_y).
 * The tilted coordinate system is the sector coordinate system rotated by 25 degrees so that the z
 * axis is perperdicular to the wire planes.
 * @author ziegler
 */
public class StateVec extends Matrix {

    private static final long serialVersionUID = 1874984192960130771L;

    private double _PathLength; // NOTE: What is this?
	private int _planeIdx;      // The wire plane index in the series of planes used in the traj.
	private double z;           // NOTE: What is this?
	private double b;           // NOTE: What is this?
	private double h;           // KF projector

    public double x() {return(get(0,0));}
    public double y() {return(get(1,0));}
    public double tanThetaX() {return(get(2,0));}
    public double tanThetaY() {return(get(3,0));}

    public double getPathLength() {return _PathLength;}
    public void setPathLength(double _PathLength) {this._PathLength = _PathLength;}

	public int getPlaneIdx() {return _planeIdx;}
    public void setPlaneIdx(int _planeIdx) {this._planeIdx = _planeIdx;}

	public double getZ() {return z;}
	public void setZ(double z) {this.z = z;}

	public double getB() {return b;}
	public void setB(double b) {this.b = b;}

	public double getProjector() {return h;}
	public void setProjector(double h) {this.h = h;}

    /** Instantiates a new, empty vector. */
    public StateVec() {
        super(4,1);
    }

    /**
     * Instantiates a new state vector with given explicit parameters.
     * @param y      the given _y
     * @param x      the given _x
     * @param tanThX the given _tanThetaX
     * @param tanThY the given _tanThetaY
     */
    public StateVec(double x, double y, double tanThX, double tanThY) {
        super(4,1);
        set(0, 0, x);
        set(1, 0, y);
        set(2, 0, tanThX);
        set(3, 0, tanThY);
    }

    /**
     * Instantiates a new state vector copying another one.
     * @param v the other vector
     */
    public StateVec(StateVec v) {
        super(4,1);
        set(0, 0, v.x());
        set(1, 0, v.y());
        set(2, 0, v.tanThetaX());
        set(3, 0, v.tanThetaY());
    }

    /**
     * Instantiates a new state vector from a Jama matrix. Needed since Jama.Matrix cannot be casted
     * into StateVec by default.
     * @param m the m
     */
    private StateVec(Matrix m) {
        super(4,1);
        set(0, 0, m.get(0, 0));
        set(1, 0, m.get(1, 0));
        set(2, 0, m.get(2, 0));
        set(3, 0, m.get(3, 0));
    }

    /**
     * Sets the stateVec with given explicit parameters.
     * @param y      the given _y
     * @param x      the given _x
     * @param tanThX the given _tanThetaX
     * @param tanThY the given _tanThetaY
     */
    public void set(double x, double y, double tanThX, double tanThY) {
        set(0, 0, x);
        set(1, 0, y);
        set(2, 0, tanThX);
        set(3, 0, tanThY);
    }

    /**
     * Sets the vector with the parameters of another vector.
     * @param V the other vector
     */
    public void set(StateVec V) {
        set(0, 0, V.x());
        set(1, 0, V.y());
        set(2, 0, V.tanThetaX());
        set(3, 0, V.tanThetaY());
    }

	/**
	 * Writes and returns a string containing the state vector's information.
	 * @return state vector's information encoded in a string.
	 */
    public String printInfo() {
		return "State Vector:" +
			   "\nx         : " + this.x() +
			   "\ny         : " + this.y() +
			   "\ntanThetaX : " + this.tanThetaX() +
			   "\ntanThetaY : " + this.tanThetaY();
    }
}
