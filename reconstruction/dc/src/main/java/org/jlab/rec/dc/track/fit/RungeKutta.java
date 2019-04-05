/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.jlab.rec.dc.track.fit;

import org.ejml.simple.SimpleMatrix;
import java.util.ArrayList;
import org.jlab.clas.swimtools.Swim;
import org.jlab.rec.dc.Constants;

/**
 *
 * @author ziegler
 */
public class RungeKutta {

    private final float[] _b = new float[3];
    final double v = 0.0029979245;
    private final ArrayList<Double> k1,  k2,  k3,  k4;
    private final ArrayList<Double> jk1, jk2, jk3, jk4;

    public RungeKutta() {
        this.k1  = new ArrayList<Double>(4);
        this.k2  = new ArrayList<Double>(4);
        this.k3  = new ArrayList<Double>(4);
        this.k4  = new ArrayList<Double>(4);
        this.jk1 = new ArrayList<Double>(12);
        this.jk2 = new ArrayList<Double>(12);
        this.jk3 = new ArrayList<Double>(12);
        this.jk4 = new ArrayList<Double>(12);
    }

    void RK4transport(int sector, double q, double x0, double y0, double z0, double tx0, double ty0,
                      double h, Swim swimmer, StateVecs.CovMat covMat, StateVecs.StateVec fVec,
                      StateVecs.CovMat fCov, double mass, double dPath) {

        // Jacobian
        double[][] u = new double[5][5];
        double[][] C = new double[5][5];
        double deltx_deltx0_0 = 1;
        double deltx_delty0_0 = 0;
        double deltx_delq0_0  = 0;
        double delty_deltx0_0 = 0;
        double delty_delty0_0 = 1;
        double delty_delq0_0  = 0;

        // Auxiliary variables
        double qv = q*v;
        double hh = 0.5*h;

        // ==- K1 -=================================================================================
        swimmer.Bfield(sector, x0, y0, z0, _b);

        // Auxiliary variables
        double C1sq = 1 + tx0*tx0 + ty0*ty0;
        double C1   = Math.sqrt(C1sq);
        double Ax1  = Ax(tx0, ty0, _b[0], _b[1], _b[2], C1);
        double Ay1  = Ay(tx0, ty0, _b[0], _b[1], _b[2], C1);

        // State
        double x1  = tx0;
        double y1  = ty0;
        double tx1 = qv * Ax1;
        double ty1 = qv * Ay1;

        // Jacobian
        double delx_deltx0_1 = deltx_deltx0_0;
        double delx_delty0_1 = deltx_delty0_0;
        double delx_delq0_1  = deltx_delq0_0;

        double dely_deltx0_1 = delty_deltx0_0;
        double dely_delty0_1 = delty_delty0_0;
        double dely_delq0_1  = delty_delq0_0;

        double deltx_deltx0_1 = qv * delAx_deltx(tx0, ty0, _b[0], _b[1], _b[2], C1sq, C1, Ax1, Ay1);
        double deltx_delty0_1 = qv * delAx_delty(tx0, ty0, _b[0], _b[1], _b[2], C1sq, C1, Ax1, Ay1);
        double deltx_delq0_1  = v  * Ax1;

        double delty_deltx0_1 = qv * delAy_deltx(tx0, ty0, _b[0], _b[1], _b[2], C1sq, C1, Ax1, Ay1);
        double delty_delty0_1 = qv * delAy_delty(tx0, ty0, _b[0], _b[1], _b[2], C1sq, C1, Ax1, Ay1);
        double delty_delq0_1  = v  * Ay1;

        // ==- K2 -=================================================================================
        swimmer.Bfield(sector, x0+hh*x1, y0+hh*y1, z0+hh, _b);

        // State 1
        double x2 = tx0 + hh*tx1;
        double y2 = ty0 + hh*ty1;

        // Auxiliary variables
        double C2sq = 1 + x2*x2 + y2*y2;
        double C2   = Math.sqrt(C2sq);
        double Ax2  = Ax(x2, y2, _b[0], _b[1], _b[2], C2);
        double Ay2  = Ay(x2, y2, _b[0], _b[1], _b[2], C2);

        double dtxtx1 = deltx_deltx0_0 + hh*deltx_deltx0_1;
        double dtxty1 = deltx_delty0_0 + hh*deltx_delty0_1;
        double dtxtq1 = deltx_delq0_0  + hh*deltx_delq0_1;
        double dtytx1 = delty_deltx0_0 + hh*delty_deltx0_1;
        double dtyty1 = delty_delty0_0 + hh*delty_delty0_1;
        double dtytq1 = delty_delq0_0  + hh*delty_delq0_1;

        // State 2
        double tx2 = qv * Ax2;
        double ty2 = qv * Ay2;

        // Jacobian
        double delx_deltx0_2 = deltx_deltx0_0 + hh*deltx_deltx0_1;
        double delx_delty0_2 = deltx_delty0_0 + hh*deltx_delty0_1;
        double delx_delq0_2  = deltx_delq0_0  + hh*deltx_delq0_1;

        double dely_deltx0_2 = delty_deltx0_0 + hh*delty_deltx0_1;
        double dely_delty0_2 = delty_delty0_0 + hh*delty_delty0_1;
        double dely_delq0_2  = delty_delq0_0  + hh*delty_delq0_1;

        double deltx_deltx0_2 = this.deltx_deltx0_next(qv, x2, y2, _b[0], _b[1], _b[2],
                                                       dtxtx1, dtytx1, C2sq, C2, Ax2, Ay2);
        double deltx_delty0_2 = this.deltx_delty0_next(qv, x2, y2, _b[0], _b[1], _b[2],
                                                       dtxty1, dtyty1, C2sq, C2, Ax2, Ay2);
        double deltx_delq0_2  = this.deltx_delq0_next(qv, v, x2, y2, _b[0], _b[1], _b[2],
                                                      dtxtq1, dtytq1, C2sq, C2, Ax2, Ay2);

        double delty_deltx0_2 = this.delty_deltx0_next(qv, x2, y2, _b[0], _b[1], _b[2],
                                                       dtxtx1, dtytx1, C2sq, C2, Ax2, Ay2);
        double delty_delty0_2 = this.delty_delty0_next(qv,x2,y2,_b[0],_b[1],_b[2],
                                                       dtxty1, dtyty1, C2sq, C2, Ax2, Ay2);
        double delty_delq0_2  = this.delty_delq0_next(qv,v,x2,y2,_b[0],_b[1],_b[2],
                                                      dtxtq1, dtytq1, C2sq, C2, Ax2, Ay2);

        // ==- K3 -=================================================================================
        swimmer.Bfield(sector, x0+hh*x2, y0+hh*y2, z0+hh, _b);

        // State 1
        double x3 = tx0 + hh*tx2;
        double y3 = ty0 + hh*ty2;

        // Auxiliary variables
        double C3sq = 1 + x3*x3 + y3*y3;
        double C3   = Math.sqrt(C3sq);
        double Ax3  = Ax(x3, y3, _b[0], _b[1], _b[2], C3);
        double Ay3  = Ay(x3, y3, _b[0], _b[1], _b[2], C3);

        double dtxtx2 = deltx_deltx0_0 + hh*deltx_deltx0_2;
        double dtxty2 = deltx_delty0_0 + hh*deltx_delty0_2;
        double dtxtq2 = deltx_delq0_0  + hh*deltx_delq0_2;
        double dtytx2 = delty_deltx0_0 + hh*delty_deltx0_2;
        double dtyty2 = delty_delty0_0 + hh*delty_delty0_2;
        double dtytq2 = delty_delq0_0  + hh*delty_delq0_2;

        // State 2
        double tx3 = qv * Ax3;
        double ty3 = qv * Ay3;

        // Jacobian
        double delx_deltx0_3 = deltx_deltx0_0 + hh*deltx_deltx0_2;
        double delx_delty0_3 = deltx_delty0_0 + hh*deltx_delty0_2;
        double delx_delq0_3  = deltx_delq0_0  + hh*deltx_delq0_2;

        double dely_deltx0_3 = delty_deltx0_0 + hh*delty_deltx0_2;
        double dely_delty0_3 = delty_delty0_0 + hh*delty_delty0_2;
        double dely_delq0_3  = delty_delq0_0  + hh*delty_delq0_2;

        double deltx_deltx0_3 = this.deltx_deltx0_next(qv, x3, y3, _b[0], _b[1], _b[2],
                                                       dtxtx2, dtytx2, C3sq, C3, Ax3, Ay3);
        double deltx_delty0_3 = this.deltx_delty0_next(qv,x3,y3,_b[0],_b[1],_b[2],
                                                       dtxty2, dtyty2, C3sq, C3, Ax3, Ay3);
        double deltx_delq0_3  = this.deltx_delq0_next(qv,v,x3,y3,_b[0],_b[1],_b[2],
                                                      dtxtq2, dtytq2, C3sq, C3, Ax3, Ay3);

        double delty_deltx0_3 = this.delty_deltx0_next(qv,x3,y3,_b[0],_b[1],_b[2],
                                                       dtxtx2, dtytx2, C3sq, C3, Ax3, Ay3);
        double delty_delty0_3 = this.delty_delty0_next(qv,x3,y3,_b[0],_b[1],_b[2],
                                                       dtxty2, dtyty2, C3sq, C3, Ax3, Ay3);
        double delty_delq0_3 = this.delty_delq0_next(qv,v,x3,y3,_b[0],_b[1],_b[2],
                                                     dtxtq2, dtytq2, C3sq, C3, Ax3, Ay3);

        // ==- K4 -=================================================================================
        swimmer.Bfield(sector, x0+h*x3, y0+h*y3, z0+h, _b);

        // State 1
        double x4 = tx0 + h*tx3;
        double y4 = ty0 + h*ty3;

        // Auxiliary variables
        double C4sq = 1 + x4*x4 + y4*y4;
        double C4   = Math.sqrt(C4sq);
        double Ax4  = Ax(x4, y4, _b[0], _b[1], _b[2], C4);
        double Ay4  = Ay(x4, y4, _b[0], _b[1], _b[2], C4);

        double dtxtx3 = deltx_deltx0_0 + hh*deltx_deltx0_3;
        double dtxty3 = deltx_delty0_0 + hh*deltx_delty0_3;
        double dtxtq3 = deltx_delq0_0  + hh*deltx_delq0_3;
        double dtytx3 = delty_deltx0_0 + hh*delty_deltx0_3;
        double dtyty3 = delty_delty0_0 + hh*delty_delty0_3;
        double dtytq3 = delty_delq0_0  + hh*delty_delq0_3;

        // State 2
        double tx4 = qv * Ax4;
        double ty4 = qv * Ay4;

        // Jacobian
        double delx_deltx0_4 = deltx_deltx0_0 + h*deltx_deltx0_3;
        double delx_delty0_4 = deltx_delty0_0 + h*deltx_delty0_3;
        double delx_delq0_4  = deltx_delq0_0  + h*deltx_delq0_3;

        double dely_deltx0_4 = delty_deltx0_0 + h*delty_deltx0_3;
        double dely_delty0_4 = delty_delty0_0 + h*delty_delty0_3;
        double dely_delq0_4  = delty_delq0_0  + h*delty_delq0_3;

        double deltx_deltx0_4 = this.deltx_deltx0_next(qv, x4, y4, _b[0], _b[1], _b[2],
                                                       dtxtx3, dtytx3, C4sq, C4, Ax4, Ay4);
        double deltx_delty0_4 = this.deltx_delty0_next(qv,x4,y4,_b[0],_b[1],_b[2],
                                                       dtxty3, dtyty3, C4sq, C4, Ax4, Ay4);
        double deltx_delq0_4  = this.deltx_delq0_next(qv,v,x4,y4,_b[0],_b[1],_b[2],
                                                      dtxtq3, dtytq3, C4sq, C4, Ax4, Ay4);

        double delty_deltx0_4 = this.delty_deltx0_next(qv,x4,y4,_b[0],_b[1],_b[2],
                                                       dtxtx3, dtytx3, C4sq, C4, Ax4, Ay4);
        double delty_delty0_4 = this.delty_delty0_next(qv,x4,y4,_b[0],_b[1],_b[2],
                                                       dtxty3, dtyty3, C4sq, C4, Ax4, Ay4);
        double delty_delq0_4  = this.delty_delq0_next(qv,v,x4,y4,_b[0],_b[1],_b[2],
                                                      dtxtq3, dtytq3, C4sq, C4, Ax4, Ay4);

        // ==- RK4 -================================================================================
        // State
        double z  = z0  + h;
        double x  = x0  + this.RK4(x1,  x2,  x3,  x4,  h);
        double y  = y0  + this.RK4(y1,  y2,  y3,  y4,  h);
        double tx = tx0 + this.RK4(tx1, tx2, tx3, tx4, h);
        double ty = ty0 + this.RK4(ty1, ty2, ty3, ty4, h);

        // Jacobian
        double delx_deltx0 = this.RK4(delx_deltx0_1, delx_deltx0_2, delx_deltx0_3, delx_deltx0_4, h);
        double delx_delty0 = this.RK4(delx_delty0_1, delx_delty0_2, delx_delty0_3, delx_delty0_4, h);
        double delx_delq0  = this.RK4(delx_delq0_1,  delx_delq0_2,  delx_delq0_3,  delx_delq0_4,  h);

        double dely_deltx0 = this.RK4(dely_deltx0_1, dely_deltx0_2, dely_deltx0_3, dely_deltx0_4, h);
        double dely_delty0 = this.RK4(dely_delty0_1, dely_delty0_2, dely_delty0_3, dely_delty0_4, h);
        double dely_delq0  = this.RK4(dely_delq0_1,  dely_delq0_2,  dely_delq0_3,  dely_delq0_4,  h);

        double deltx_deltx0 = this.RK4(deltx_deltx0_1, deltx_deltx0_2, deltx_deltx0_3, deltx_deltx0_4, h) + 1;
        double deltx_delty0 = this.RK4(deltx_delty0_1, deltx_delty0_2, deltx_delty0_3, deltx_delty0_4, h);
        double deltx_delq0  = this.RK4(deltx_delq0_1,  deltx_delq0_2,  deltx_delq0_3,  deltx_delq0_4,  h);

        double delty_deltx0 = this.RK4(delty_deltx0_1, delty_deltx0_2, delty_deltx0_3, delty_deltx0_4, h);
        double delty_delty0 = this.RK4(delty_delty0_1, delty_delty0_2, delty_delty0_3, delty_delty0_4, h) + 1;
        double delty_delq0  = this.RK4(delty_delq0_1,  delty_delq0_2,  delty_delq0_3,  delty_delq0_4,  h);

        // ==- Covariance Matrix -==================================================================
        // covMat = FCF^T; u = FC;
        for (int j1 = 0; j1 < 5; j1++) {
            u[0][j1] = covMat.covMat.get(0,j1) + covMat.covMat.get(2,j1) * delx_deltx0  + covMat.covMat.get(3,j1) * delx_delty0 + covMat.covMat.get(4,j1) * delx_delq0;
            u[1][j1] = covMat.covMat.get(1,j1) + covMat.covMat.get(2,j1) * dely_deltx0  + covMat.covMat.get(3,j1) * dely_delty0 + covMat.covMat.get(4,j1) * dely_delq0;
            u[2][j1] = covMat.covMat.get(2,j1) + covMat.covMat.get(3,j1) * deltx_delty0 + covMat.covMat.get(4,j1) * deltx_delq0;
            u[3][j1] = covMat.covMat.get(2,j1) * delty_deltx0 + covMat.covMat.get(3,j1) + covMat.covMat.get(4,j1) * delty_delq0;
            u[4][j1] = covMat.covMat.get(4,j1);
        }

        for (int i1 = 0; i1 < 5; i1++) {
            C[i1][0] = u[i1][0] + u[i1][2] * delx_deltx0 + u[i1][3] * delx_delty0 + u[i1][4] * delx_delq0;
            C[i1][1] = u[i1][1] + u[i1][2] * dely_deltx0 + u[i1][3] * dely_delty0 + u[i1][4] * dely_delq0;
            C[i1][2] = u[i1][2] + u[i1][3] * deltx_delty0 + u[i1][4] * deltx_delq0;
            C[i1][3] = u[i1][2] * delty_deltx0 + u[i1][3] + u[i1][4] * delty_delq0;
            C[i1][4] = u[i1][4];
        }

        // ==- Q Process Noise Matrix Estimate -====================================================
        double p  = Math.abs(1. / q);
        double pz = p / Math.sqrt(1 + tx * tx + ty * ty);
        double px = tx * pz;
        double py = ty * pz;

        // Path length in radiation length units = t/X0 [true path length/ X0]
        // Ar radiation length = 14 cm
        double t_ov_X0 = Math.signum(h) * h / Constants.ARGONRADLEN;

        // Use particle momentum
        double beta = p / Math.sqrt(p * p + mass * mass);
        double cosEntranceAngle = Math.abs((x * px + y * py + z * pz) / (Math.sqrt(x * x + y * y + z * z) * p));
        double pathLength = t_ov_X0 / cosEntranceAngle;

        // Highland-Lynch-Dahl formula
        double sctRMS = (0.0136 / (beta * p)) * Math.sqrt(pathLength) * (1 + 0.038 * Math.log(pathLength));

        double cov_txtx = (1 + tx * tx) * (1 + tx * tx + ty * ty) * sctRMS * sctRMS;
        double cov_tyty = (1 + ty * ty) * (1 + tx * tx + ty * ty) * sctRMS * sctRMS;
        double cov_txty = tx * ty * (1 + tx * tx + ty * ty) * sctRMS * sctRMS;

        if (h > 0) {
            C[2][2] += cov_txtx;
            C[2][3] += cov_txty;
            C[3][2] += cov_txty;
            C[3][3] += cov_tyty;
        }

        fVec.x  = x;
        fVec.y  = y ;
        fVec.z  = z0+h;
        fVec.tx = tx;
        fVec.ty = ty;
        fVec.Q  = q;
        fVec.B  = Math.sqrt(_b[0]*_b[0]+_b[1]*_b[1]+_b[2]*_b[2]);
        fVec.deltaPath = Math.sqrt((x0-x)*(x0-x)+(y0-y)*(y0-y)+h*h)+dPath;
        fCov.covMat    = new SimpleMatrix(C);
    }

    private double RK4(double k1, double k2, double k3, double k4, double h) {
        return h/6 * (k1 + 2*k2 + 2*k3 + k4);
    }

    private double Ax(double tx, double ty, double Bx, double By, double Bz, double C) {
        return C * (ty * (tx * Bx + Bz) - (1 + tx * tx) * By);
    }
    private double Ay(double tx, double ty, double Bx, double By, double Bz, double C) {
        return C * (-tx * (ty * By + Bz) + (1 + ty * ty) * Bx);
    }

    private double delAx_deltx(double tx, double ty, double Bx, double By, double Bz,
                               double Csq, double C, double Ax, double Ay) {
        return tx * Ax / Csq + C * (ty * Bx - 2 * tx * By);
    }
    private double delAx_delty(double tx, double ty, double Bx, double By, double Bz,
                               double Csq, double C, double Ax, double Ay) {
        return ty * Ax / Csq + C * (tx * Bx + Bz);
    }

    private double delAy_deltx(double tx, double ty, double Bx, double By, double Bz,
                               double Csq, double C, double Ax, double Ay) {
        return tx * Ay / Csq + C * (-ty * By - Bz);
    }
    private double delAy_delty(double tx, double ty, double Bx, double By, double Bz,
                               double Csq, double C, double Ax, double Ay) {
        return ty * Ay / Csq + C * (-tx * By + 2 * ty * Bx);
    }

    private double deltx_deltx0_next(double qv, double tx1, double ty1, float b0, float b1, float b2,
                                     double deltx_deltx0_1, double delty_deltx0_1,
                                     double Csq, double C, double Ax, double Ay) {
        return qv * (delAx_deltx(tx1, ty1, b0, b1, b2, Csq, C, Ax, Ay) * deltx_deltx0_1
                     + delAx_delty(tx1, ty1, b0, b1, b2, Csq, C, Ax, Ay) * delty_deltx0_1);
    }
    private double deltx_delty0_next(double qv, double tx1, double ty1, float b0, float b1, float b2,
                                     double deltx_delty0_1, double delty_delty0_1,
                                     double Csq, double C, double Ax, double Ay) {
        return qv * (delAx_deltx(tx1, ty1, b0, b1, b2, Csq, C, Ax, Ay) * deltx_delty0_1 // NOTE: HERE
                     + delAx_delty(tx1, ty1, b0, b1, b2, Csq, C, Ax, Ay) * delty_delty0_1);
    }
    private double deltx_delq0_next(double qv, double v, double tx1, double ty1, float b0, float b1, float b2,
                                    double deltx_delq0_1, double delty_delq0_1,
                                    double Csq, double C, double Ax, double Ay) {
        return v * Ax + qv * (delAx_deltx(tx1, ty1, b0, b1, b2, Csq, C, Ax, Ay) * deltx_delq0_1
                              + delAx_delty(tx1, ty1, b0, b1, b2, Csq, C, Ax, Ay) * delty_delq0_1);
    }

    private double delty_deltx0_next(double qv, double tx1, double ty1, float b0, float b1, float b2,
                                     double deltx_deltx0_1, double delty_deltx0_1,
                                     double Csq, double C, double Ax, double Ay) {
        return qv * (delAy_deltx(tx1, ty1, b0, b1, b2, Csq, C, Ax, Ay) * deltx_deltx0_1
                     + delAy_delty(tx1, ty1, b0, b1, b2, Csq, C, Ax, Ay) * delty_deltx0_1);
    }
    private double delty_delty0_next(double qv, double tx1, double ty1, float b0, float b1, float b2,
                                     double deltx_delty0_1, double delty_delty0_1,
                                     double Csq, double C, double Ax, double Ay) {
        return qv * (delAy_deltx(tx1, ty1, b0, b1, b2, Csq, C, Ax, Ay) * deltx_delty0_1
                     + delAy_delty(tx1, ty1, b0, b1, b2, Csq, C, Ax, Ay) * delty_delty0_1);
    }
    private double delty_delq0_next(double qv, double v, double tx1, double ty1, float b0, float b1, float b2,
                                    double deltx_delq0_1, double delty_delq0_1,
                                    double Csq, double C, double Ax, double Ay) {
        return v * Ay + qv * (delAy_deltx(tx1, ty1, b0, b1, b2, Csq, C, Ax, Ay) * deltx_delq0_1
                              + delAy_delty(tx1, ty1, b0, b1, b2, Csq, C, Ax, Ay) * delty_delq0_1);
    }
}
