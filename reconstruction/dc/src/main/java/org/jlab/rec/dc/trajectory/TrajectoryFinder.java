package org.jlab.rec.dc.trajectory;

import java.util.ArrayList;
import java.util.List;

import org.jlab.clas.swimtools.Swim;
import org.jlab.geom.prim.Vector3D;
import trackfitter.fitter.LineFitter;
import org.jlab.detector.geant4.v2.DCGeant4Factory;

import org.jlab.rec.dc.Constants;
import org.jlab.rec.dc.cross.Cross;

/**
 * A driver class to find the trajectory of a track candidate. THE PATH TO FIELD MAPS IS SET BY THE
 * CLARA_SERVICES ENVIRONMENT VAR.
 * @author ziegler
 */
public class TrajectoryFinder {

	private LineFitter lineFit;

	private double[] tanTheta_x_fitCoeff;
	private double[] tanTheta_y_fitCoeff;
	private double[] x_fitCoeff;
	private double[] y_fitCoeff;

	/** Step size used in integral Bdl Riemann integration. */
	public double mmStepSizeForIntBdl = 10;

	private double PathLength;

	int counter = 0;
    public double TrajChisqProbFitXZ;
    public double TrajChisqProbFitYZ;

    /**
     * Determines and returns a trajectory from an input list of crosses.
     * @param candCrossList the list of crosses
     * @return              a trajectory object
     */
    public Trajectory findTrajectory(List<Cross> candCrossList,
									 DCGeant4Factory DcDetector,
									 Swim dcSwim) {

		Trajectory traj = new Trajectory();
        if (candCrossList.isEmpty()) return traj;

        traj.addAll(candCrossList);
        traj.set_Sector(candCrossList.get(0).get_Sector());
        fitTrajectory(traj);
        if (this.TrajChisqProbFitXZ < Constants.TCHISQPROBFITXZ) return null;

        traj.set_Trajectory(getStateVecsAlongTrajectory(DcDetector));
        traj.set_IntegralBdl(integralBdl(candCrossList.get(0).get_Sector(), DcDetector, dcSwim));
        traj.set_PathLength(PathLength);

        return traj;
    }

    /**
     * NOTE: Lacks JavaDoc description.
     * @param DcDetector DC detector utility
     * @return           integral Bdl
     */
    public double integralBdl(int sector, DCGeant4Factory DcDetector, Swim dcSwim) {

        double z1 = DcDetector.getRegionMidpoint(0).z;
        double z3 = DcDetector.getRegionMidpoint(2).z;

        double z  = z1;

        double intBdl  = 0;
        double pathLen = 0;
        double x0 = x_fitCoeff[0] * z1*z1 + x_fitCoeff[1] * z1 + x_fitCoeff[2];
        double y0 = y_fitCoeff[0] * z1*z1 + y_fitCoeff[1] * z1 + y_fitCoeff[2];
        double z0 = z1;


        while(z <= z3) {
            counter++;
            double x = x_fitCoeff[0] * z*z + x_fitCoeff[1] * z + x_fitCoeff[2];
            double y = y_fitCoeff[0] * z*z + y_fitCoeff[1] * z + y_fitCoeff[2];

            float[] result = new float[3];
            dcSwim.Bfield(sector, (x + x0) * 0.5, (y + y0) * 0.5, (z + z0) * 0.5, result);

			Vector3D dl = new Vector3D(x - x0, 0, z - z0);
            Vector3D Bf = new Vector3D(result[0], result[1], result[2]);

			intBdl  += dl.cross(Bf).mag();
            pathLen += dl.mag();
			x0 = x;
            y0 = y;
            z0 = z;

            z += mmStepSizeForIntBdl / 10.;
        }
        PathLength = pathLen;

        return intBdl;
    }

    /**
     * Determines and returns a list of state vectors along the trajectory.
     * @return the list of state vectors
     */
    public List<StateVec> getStateVecsAlongTrajectory(DCGeant4Factory DcDetector) {

		List<StateVec> stateVecAtPlanesList = new ArrayList<StateVec>(36);
        for (int superlayerIdx = 0; superlayerIdx < 6; superlayerIdx++) {
            for (int layerIdx = 0; layerIdx < 6; layerIdx++) {

				double z = DcDetector.getLayerMidpoint(superlayerIdx, layerIdx).z;
                double x = x_fitCoeff[0] * z*z + x_fitCoeff[1] * z + x_fitCoeff[2];
                double y = y_fitCoeff[0] * z*z + y_fitCoeff[1] * z + y_fitCoeff[2];
                double tanTheta_x = x_fitCoeff[0] * z + x_fitCoeff[1];
                double tanTheta_y = y_fitCoeff[0] * z + y_fitCoeff[1];

                StateVec stateVec = new StateVec(x,y,tanTheta_x, tanTheta_y);
                stateVecAtPlanesList.add(stateVec);
            }
        }
        return stateVecAtPlanesList;
    }

    /**
     * Gets the parametric form of the trajectory determined from fitting the tangent values of the
	 * state vectors linearly and constraining the quadratic parameters of the function describing
	 * the position values of these vectors.
     * @param candCrossList list of crosses used in the fit
     */
    public void fitTrajectory(List<Cross> candCrossList) {

		tanTheta_x_fitCoeff = new double[2];
        tanTheta_y_fitCoeff = new double[2];
        x_fitCoeff          = new double[3];
        y_fitCoeff          = new double[3];

        double[] theta_x     = new double[3];
        double[] theta_x_err = new double[3];
        double[] theta_y     = new double[3];
        double[] theta_y_err = new double[3];

        double[] x     = new double[3];
        double[] x_err = new double[3];
        double[] y     = new double[3];
        double[] y_err = new double[3];
        double[] z     = new double[3];

        for (int i = 0; i < 3; i++) {
            // Make sure that the track direction makes sense
            if (candCrossList.get(i).get_Dir().z() == 0) return;

            x[i]     = candCrossList.get(i).get_Point().x();
            x_err[i] = candCrossList.get(i).get_PointErr().x();
            y[i]     = candCrossList.get(i).get_Point().y();
            y_err[i] = candCrossList.get(i).get_PointErr().y();
            z[i]     = candCrossList.get(i).get_Point().z();

            theta_x[i] = candCrossList.get(i).get_Dir().x() / candCrossList.get(i).get_Dir().z();
            theta_x_err[i] = calcTanErr(candCrossList.get(i).get_Dir().x(),
										candCrossList.get(i).get_Dir().z(),
										candCrossList.get(i).get_DirErr().x(),
										candCrossList.get(i).get_DirErr().z());

			theta_y[i] = candCrossList.get(i).get_Dir().y() / candCrossList.get(i).get_Dir().z();
            theta_y_err[i] = calcTanErr(candCrossList.get(i).get_Dir().y(),
										candCrossList.get(i).get_Dir().z(),
										candCrossList.get(i).get_DirErr().y(),
										candCrossList.get(i).get_DirErr().z());
        }

        lineFit = new LineFitter();
        boolean linefitstatusOK = lineFit.fitStatus(z, theta_x, new double[3], theta_x_err, 3);

        // tan_thetax = alpha*z + beta;
        // x = a*z^2 +b*z +c
        if (linefitstatusOK) {
            double alpha = lineFit.getFit().slope();
            double beta  = lineFit.getFit().intercept();

            double a = alpha/2;
            double b = beta;

            double sum_inv_xerr  = 0;
            double sum_X_ov_errX = 0;

            for (int i = 0; i < 3; i++) {
                x[i] -= a * z[i]*z[i] + b*z[i];
                sum_inv_xerr  += 1. / x_err[i];
                sum_X_ov_errX += x[i] / x_err[i];
            }
            if (sum_inv_xerr == 0) return;

            double c = sum_X_ov_errX / sum_inv_xerr;

            tanTheta_x_fitCoeff[0] = alpha;
            tanTheta_x_fitCoeff[1] = beta;

            x_fitCoeff[0] = a;
            x_fitCoeff[1] = b;
            x_fitCoeff[2] = c;
        }

        TrajChisqProbFitXZ = lineFit.getFit().getProb();

        lineFit = new LineFitter();
        linefitstatusOK = lineFit.fitStatus(z, theta_y, new double[3], theta_y_err, 3);

        // tan_thetay = alpha*z + beta;
        // y = a*z^2 +b*z + c
        if (linefitstatusOK) {
            double alpha = lineFit.getFit().slope();
            double beta = lineFit.getFit().intercept();

            double a = alpha/2;
            double b = beta;

            double sum_inv_yerr  = 0;
            double sum_Y_ov_errY = 0;
            for (int i = 0; i < 3; i++) {
                y[i] -= a * z[i]*z[i] + b*z[i];
                sum_inv_yerr  += 1. / y_err[i];
                sum_Y_ov_errY += y[i] / y_err[i];
            }

			if (sum_inv_yerr == 0) return;

            double c = sum_Y_ov_errY/sum_inv_yerr;

            tanTheta_y_fitCoeff[0] = alpha;
            tanTheta_y_fitCoeff[1] = beta;

            y_fitCoeff[0] = a;
            y_fitCoeff[1] = b;
            y_fitCoeff[2] = c;
        }

        TrajChisqProbFitYZ = lineFit.getFit().getProb();
    }

    /**
     * Calculates the error on the tangent num/denom.
     * @param num     NOTE: Missing description
     * @param denom   NOTE: Missing description
     * @param numEr   NOTE: Missing description
     * @param denomEr NOTE: Missing description
     * @return        the error
     */
    private double calcTanErr(double num, double denom, double numEr, double denomEr) {
        double d1 = num / (denom*denom);
        double e1 = denomEr;
        double d2 = 1 / denom;
        double e2 = numEr;

        return Math.sqrt(d1*d1 + e1*e1 + d2*d2 * e2*e2);
    }
}
