package org.jlab.rec.dc.timetodistance;

import java.math.RoundingMode;
import java.text.DecimalFormat;

import org.jlab.detector.calib.utils.DatabaseConstantProvider;
import org.jlab.utils.groups.IndexedTable;

import org.jlab.rec.dc.Constants;

/**
 * @author ziegler
 */
public class TableLoader {

    public TableLoader() {}

    static final protected int nBinsT = 2000;
    public static double[][][][][] DISTFROMTIME
            = new double[  6   ][ 6  ][  8  ][  6   ][nBinsT];
    //                  [sector][slyr][alpha][Bfield][time][bins]
    static boolean T2DLOADED  = false;
    static boolean T0LOADED   = false;
    static int minBinIdxB     = 0;
    static int maxBinIdxB     = 7;
    static int minBinIdxAlpha = 0;
    static int maxBinIdxAlpha = 6;
    static int minBinIdxT     = 0;
    static int[][][][] maxBinIdxT  = new int[6][6][8][6];

    // Fraction of dmax corresponding to the point in the cell where the velocity is minimal
    public static double FracDmaxAtMinVel = 0.615;

    public void test(){
        TimeToDistanceEstimator tde = new TimeToDistanceEstimator();
        for (int s = 0; s < 1; s++) { // Loop over sectors
            for (int r = 2; r < 3; r++) { // Loop over slys
                for (int ibfield = 0; ibfield < 8; ibfield++) {
                    for (int icosalpha = 0; icosalpha < 6; icosalpha++) {
                        for (int tb = 0; tb < maxBinIdxT[s][r][ibfield][icosalpha]; tb++) {
                            // Math.cos(Math.toRadians(30.)) = 0.8660254037844387
                            // (1. - Math.cos(Math.toRadians(30.)))/5. = 0.02679491924311226
                            double Xalpha = -(Math.toDegrees(Math.acos(0.8660254037844387
                                    + (icosalpha) * 0.02679491924311226)) - 30.);
                            double Xtime = (2*tb + 1);
                            double Xdoca = tde.interpolateOnGrid((double) ibfield * 0.5,
                                                                 Xalpha, Xtime, s, r);

                            System.out.println(
                                  "s         : " + (s+1) +
                                "\nsl        : " + (r+1) +
                                "\ntime      : " + (2*tb+1) +
                                "\nicosalpha : " + icosalpha +
                                "\nXalpha    : " + Xalpha +
                                "\nB         : " + ibfield*0.5 +
                                "\ndis       : " + (float)DISTFROMTIME[s][r][ibfield][icosalpha][tb] +
                                "\n          : " + (float) Xdoca);
                        }
                    }
                }
            }
        }
    }

    public static synchronized void FillT0Tables(int run, String variation) {
        if (T0LOADED) return;

        System.out.println("T0 TABLE FILLED..... for Run " + run + " with VARIATION " + variation);
        DatabaseConstantProvider dbprovider = new DatabaseConstantProvider(run, variation);
        dbprovider.loadTable("/calibration/dc/time_corrections/T0Corrections");

        // Disconnect from database, it's important to do this after loading tables.
        dbprovider.disconnect();

        // T0-subtraction
        double[][][][] T0 ;
        double[][][][] T0ERR ;

        // T0s
        //                [nSec][nSL][nSlots][nCables]
        T0    = new double[6   ][6  ][7     ][6      ];
        T0ERR = new double[6   ][6  ][7     ][6      ];
        int length = dbprovider.length("/calibration/dc/time_corrections/T0Corrections/Sector");
        for (int i = 0; i < length; i++) {
            int iSec = dbprovider.getInteger("/calibration/dc/time_corrections/T0Corrections/Sector", i);
            int iSly = dbprovider.getInteger("/calibration/dc/time_corrections/T0Corrections/Superlayer", i);
            int iSlot = dbprovider.getInteger("/calibration/dc/time_corrections/T0Corrections/Slot", i);
            int iCab = dbprovider.getInteger("/calibration/dc/time_corrections/T0Corrections/Cable", i);
            double t0 = dbprovider.getDouble("/calibration/dc/time_corrections/T0Corrections/T0Correction", i);
            double t0Error = dbprovider.getDouble("/calibration/dc/time_corrections/T0Corrections/T0Error", i);

            T0[iSec - 1][iSly - 1][iSlot - 1][iCab - 1]    = t0;
            T0ERR[iSec - 1][iSly - 1][iSlot - 1][iCab - 1] = t0Error;
            Constants.setT0(T0);
            Constants.setT0Err(T0ERR);
        }

        T0LOADED = true;
    }

    public static synchronized void Fill(IndexedTable tab) {
        if (T2DLOADED) return;
        System.out.println(" T2D TABLE FILLED.....");
        double stepSize = 0.0010;
        DecimalFormat df = new DecimalFormat("#");
        df.setRoundingMode(RoundingMode.CEILING);

        for(int s = 0; s < 6; s++){ // Loop over sectors
            for(int r = 0; r < 6; r++){ // Loop over slys
                double dmax = 2.*Constants.wpdist[r];
                double tmax = tab.getDoubleValue("tmax", s + 1, r + 1,0);

                for (int ibfield = 0; ibfield < 8; ibfield++) {
                    double bfield = (double) ibfield * 0.5;
                    double maxdist = 0;

                    for (int icosalpha = 0; icosalpha < 6; icosalpha++) {
                        // Math.cos(Math.toRadians(30.)) = 0.8660254037844387
                        // 1. - Math.cos(Math.toRadians(30.)))/5. = 0.02679491924311226
                        double cos30minusalpha = 0.8660254037844387 +
                                                 (double) (icosalpha)*0.02679491924311226;

                        double alpha = -(Math.toDegrees(Math.acos(cos30minusalpha)) - 30);

                        int nxmax = (int) (dmax/stepSize);

                        for (int idist = 0; idist < nxmax; idist++) {

                            double x = (double)(idist + 1) * stepSize;
                            double timebfield = calc_Time(x, dmax, tmax, alpha, bfield, s, r, tab);

                            if (timebfield <= calc_Time(dmax, dmax, tmax, 30, bfield, s, r, tab)) {
                                maxdist = x;
                            }
                            else {
                                x = maxdist;
                            }
                            int tbin = Integer.parseInt(df.format(timebfield/2.)) -1;

                            if (tbin < 0)       tbin = 0;
                            if (tbin >= nBinsT) tbin = nBinsT - 1;
                            if (tbin > maxBinIdxT[s][r][ibfield][icosalpha]) {
                                maxBinIdxT[s][r][ibfield][icosalpha] = tbin;
                            }
                            if (DISTFROMTIME[s][r][ibfield][icosalpha][tbin] == 0) {
                                DISTFROMTIME[s][r][ibfield][icosalpha][tbin] = x;
                            } else {
                                DISTFROMTIME[s][r][ibfield][icosalpha][tbin] += stepSize;
                            }
                        }
                    }
                }
            }
        }
        T2DLOADED = true;
    }

    /**
     * Calculates the time (ns) when given inputs of distance (cm), local angle alpha (degrees) and
     * magnitude of bfield (Tesla).
     * @param x      distance to wire in cm
     * @param dmax   max distance to wire in cm
     * @param tmax   max drift time in ns
     * @param alpha  local angle in deg
     * @param bfield B field value a x in T
     * @param s      sector idx
     * @param r      superlayer idx
     * @param tab
     * @return time (ns)
     */
    public static synchronized double calc_Time(double x, double dmax, double tmax, double alpha,
            double bfield, int s, int r, IndexedTable tab) {

        // Assume a functional form (time=x/v0+a*(x/dmax)**n+b*(x/dmax)**m) for time as a function
        // of x for theta = 30 deg.

        // First, calculate n
        double deltanm = tab.getDoubleValue("deltanm", s + 1, r + 1, 0);
        double n = (1. + (deltanm - 1.) * Math.pow(FracDmaxAtMinVel, deltanm))
                 / (1. - Math.pow(FracDmaxAtMinVel, deltanm));

        // Now, calculate m
        double m = n + deltanm;

        // Determine b from the requirement that the time = tmax at dist = dmax
        double v0 = tab.getDoubleValue("v0", s + 1, r + 1, 0);
        double b  = (tmax - dmax/v0) / (1. - m/n);
        // Determine a from the requirement that the derivative at d = dmax equal the derivative
        // at d=0
        double a = -b*m/n;

        double cos30minusalpha = Math.cos(Math.toRadians(30. - alpha));
        double xhat      = x / dmax;
        double dmaxalpha = dmax * cos30minusalpha;
        double xhatalpha = x / dmaxalpha;

        // Now calculate the dist to time function for theta = 'alpha' deg.
        // Assume a functional form with the SAME POWERS N and M and coefficient a but a new
        // coefficient 'balpha' to replace b. Calculate balpha from the constraint that the value
        // of the function at dmax*cos30minusalpha is equal to tmax.

        // parameter balpha (function of the 30 degree paramters a,n,m)
        double balpha = (tmax - dmaxalpha/v0 - a * Math.pow(cos30minusalpha,n))
                      / Math.pow(cos30minusalpha, m);

        // Now calculate function
        double time = x/v0 + a*Math.pow(xhat, n) + balpha*Math.pow(xhat, m);

        // And here's a parameterization of the change in time due to a non-zero bfield for where
        // xhat = x/dmaxalpha where dmaxalpha is the 'dmax' for a track with local angle alpha (for
        // local angle = alpha)
        double delBf = tab.getDoubleValue("delta_bfield_coefficient", s + 1, r + 1, 0);

        double deltatime_bfield = delBf*Math.pow(bfield,2) * tmax *
                                (tab.getDoubleValue("b1", s + 1, r + 1, 0) * xhatalpha
                               + tab.getDoubleValue("b2", s + 1, r + 1, 0) * Math.pow(xhatalpha, 2)
                               + tab.getDoubleValue("b3", s + 1, r + 1, 0) * Math.pow(xhatalpha, 3)
                               + tab.getDoubleValue("b4", s + 1, r + 1, 0) * Math.pow(xhatalpha, 4));

        // Calculate the time at alpha degree and at a non-zero bfield
        time += deltatime_bfield;

        // Added deta(T0) correction
        time += tab.getDoubleValue("delta_T0", s + 1, r + 1, 0);

        return time;
    }
}
