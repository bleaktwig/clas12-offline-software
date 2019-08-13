package org.jlab.rec.dc.track;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

import org.jlab.clas.swimtools.Swim;

public class BFieldInterpolator {
    private Swim dcSwim;
    private int sector;
    private double[] X, Y, Z, ss;
    private HashMap<Integer, double[]> meas;

    public BFieldInterpolator(Swim dcSwim, int sector, double[] min, double[] max, double[] ss) {
        if (((max[0]-min[0])/ss[0])%1 != 0) return;
        if (((max[1]-min[1])/ss[1])%1 != 0) return;
        if (((max[2]-min[2])/ss[2])%1 != 0) return;

        this.dcSwim = dcSwim;
        this.ss     = ss;
        this.sector = sector;

        int XSize = (int) (1 + (max[0] - min[0]) / ss[0]);
        int YSize = (int) (1 + (max[1] - min[1]) / ss[1]);
        int ZSize = (int) (1 + (max[2] - min[2]) / ss[2]);

        X = new double[XSize];
        Y = new double[YSize];
        Z = new double[ZSize];

        for (int i = 0; i < XSize; ++i) X[i] = min[0] + ss[0]*i;
        for (int i = 0; i < YSize; ++i) Y[i] = min[1] + ss[1]*i;
        for (int i = 0; i < ZSize; ++i) Z[i] = min[2] + ss[2]*i;

        meas = new HashMap<Integer, double[]>();

        for (double x : X) {
            for (double y : Y) {
                for (double z : Z) {
                    float[] b = new float[3];
                    double[] db = new double[3];
                    dcSwim.Bfield(sector, x, y, z, b);

                    for (int i = 0; i < 3; ++i) db[i] = (double) b[i];
                    meas.put(this.formatInput(x,y,z), db);
                }
            }
        }
    }

    public double[] getB(double x, double y, double z) {
        if (x < X[0] || x > X[X.length-1]) return this.getRealB(x, y, z);
        if (y < Y[0] || y > Y[Y.length-1]) return this.getRealB(x, y, z);
        if (z < Z[0] || z > Z[Z.length-1]) return this.getRealB(x, y, z);

        int x0i = x == X[X.length-1] ? X.length-2 : (int) Math.floor((x - X[0])/ss[0]);
        int y0i = y == Y[Y.length-1] ? Y.length-2 : (int) Math.floor((y - Y[0])/ss[1]);
        int z0i = z == Z[Z.length-1] ? Z.length-2 : (int) Math.floor((z - Z[0])/ss[0]);

        double xd = (x - X[x0i])/ss[0];
        double yd = (y - Y[y0i])/ss[1];
        double zd = (z - Z[z0i])/ss[2];

        double[] c000 = this.getSampledB(X[x0i],   Y[y0i],   Z[z0i]);
        double[] c001 = this.getSampledB(X[x0i],   Y[y0i],   Z[z0i+1]);
        double[] c010 = this.getSampledB(X[x0i],   Y[y0i+1], Z[z0i]);
        double[] c011 = this.getSampledB(X[x0i],   Y[y0i+1], Z[z0i+1]);
        double[] c100 = this.getSampledB(X[x0i+1], Y[y0i],   Z[z0i]);
        double[] c101 = this.getSampledB(X[x0i+1], Y[y0i],   Z[z0i+1]);
        double[] c110 = this.getSampledB(X[x0i+1], Y[y0i+1], Z[z0i]);
        double[] c111 = this.getSampledB(X[x0i+1], Y[y0i+1], Z[z0i+1]);

        double[] c00 = new double[3];
        double[] c01 = new double[3];
        double[] c10 = new double[3];
        double[] c11 = new double[3];
        for (int i = 0; i < 3; ++i) {
            c00[i] = c000[i]*(1 - xd) + c100[i]*xd;
            c01[i] = c001[i]*(1 - xd) + c101[i]*xd;
            c10[i] = c010[i]*(1 - xd) + c110[i]*xd;
            c11[i] = c011[i]*(1 - xd) + c111[i]*xd;
        }

        double[] c0 = new double[3];
        double[] c1 = new double[3];
        for (int i = 0; i < 3; ++i) {
            c0[i] = c00[i]*(1 - yd) + c10[i]*yd;
            c1[i] = c01[i]*(1 - yd) + c11[i]*yd;
        }

        double[] c = new double[3];
        for (int i = 0; i < 3; ++i) c[i] = c0[i]*(1 - zd) + c1[i]*zd;

        return c;
    }

    private int formatInput(double x, double y, double z) {
        return (int)(z*1.e6 + (x-X[0])*1.e3 + (y-Y[0]));
    }

    private double[] getSampledB(double x, double y, double z) {
        return meas.get(formatInput(x, y, z));
    }

    private double[] getRealB(double x, double y, double z) {
        float[] b   = new float[3];
        // System.out.printf("🐦 - accessed B outside range! (%3.0f %3.0f %3.0f)\n", x, y, z);
        dcSwim.Bfield(sector, x, y, z, b);
        return new double[] {(double) b[0], (double) b[1], (double) b[2]};
    }
}
