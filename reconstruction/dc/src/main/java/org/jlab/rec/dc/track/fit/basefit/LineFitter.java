package org.jlab.rec.dc.track.fit.basefit;

import java.util.List;
import java.util.ArrayList;

/**
 * A least square fitting method.
 * For a linear fit, f(a,b) = a + bx taking y errors into account.
 *
 * @author ziegler
 */
public class LineFitter {

	private LineFitPars _linefitresult;

	// fit status
	public boolean fitStatus(List<Double> x,       List<Double> y,
							 List<Double> sigma_x, List<Double> sigma_y,
							 int nbpoints) {
		// must have enough points to do the fit
		if (nbpoints < 2) {
			return false;
		}

		// now do the fit
		// initialize weight-sum and moments
		double Sw, Sx, Sy, Sxx, Sxy;
		Sw = Sx = Sy = Sxx = Sxy = 0.;

		ArrayList<Double> w = new ArrayList<Double>(nbpoints);

		for (int i = 0; i < nbpoints; i++) {
			double sigmaPreCalc = (sigma_y.get(i) * sigma_y.get(i) +
						  		   sigma_x.get(i) * sigma_x.get(i));
			if (sigmaPreCalc == 0) {
				return false;
			}
			w.add(i, 1./sigmaPreCalc);
			Sw  += w.get(i);
			// the moments
			double xPreCalc = w.get(i) * x.get(i);
			Sx  += xPreCalc;
			Sy  += w.get(i) * y.get(i);
			Sxy += xPreCalc * y.get(i);
			Sxx += xPreCalc * x.get(i);
		}

		// the determinant; must be > 0
		double determ = Sw*Sxx - Sx*Sx;
		if(determ < 1e-19) determ = 1e-19; //straight track approximation

		double slopeSol  = (Sw*Sxy - Sx*Sy)/determ;
		double intercSol = (Sy*Sxx - Sx*Sxy)/determ;
		// the errors on these parameters
		double slopeErr  = Math.sqrt(Sw/determ);
		double intercErr = Math.sqrt(Sxx/determ);
		double SlInCov   = -Sx/determ;

		if (Math.abs(slopeSol) < 0 || Math.abs(intercSol) < 0) {
			return false;
		}

		// calculate the chi^2
		double chi_2 = 0.;
		// individual chi^2 for each fitted point
		double pointchi_2[] = new double[nbpoints];
		for (int i = 0; i < nbpoints; i++) {
			pointchi_2[i] = w.get(i) *
							((y.get(i) - (slopeSol * x.get(i) + intercSol)) *
							 (y.get(i) - (slopeSol * x.get(i) + intercSol)));
			chi_2 += pointchi_2[i];
		}

		// the number of degrees of freedom
		int ndf = nbpoints - 2;

		// instantiate fit object to be returned;
		_linefitresult = new LineFitPars(slopeSol,   intercSol,
										 slopeErr,   intercErr,
										 SlInCov,    chi_2,
										 pointchi_2, ndf);

		// if there is a fit return true
		return true;
	}

	// return the fit result
	public LineFitPars getFit() {
		return _linefitresult;
	}
}
