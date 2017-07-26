#include <stdio.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_errno.h>
#include "susceptibility_tensor.h"

double MJ(struct params * params)
{
	double ans = exp(-params->gamma/params->theta_e) 
		   / (4. * M_PI * params->theta_e*params->theta_e 
	  	      * gsl_sf_bessel_Kn(2, 1./params->theta_e));
	return ans;
}

double PL(struct params * params)
{
	double beta = sqrt(1. - 1./pow(params->gamma, 2.));

	double ans = (params->pl_p - 1.) * (-1 + 2. * params->gamma * params->gamma 
					 + params->pl_p * (params->gamma*params->gamma - 1.))
		    / (4. * M_PI * (pow(params->gamma_min, -1. - params->pl_p) - pow(params->gamma_max, -1. - params->pl_p))
			* beta * (params->gamma*params->gamma - 1.)) * pow(params->gamma, -3. - params->pl_p);
	return ans;	
}

double Df(struct params * params)
{
	if(params->dist == 0)
	{
		return MJ(params);
	}
	else if(params->dist == 1)
	{
		return PL(params);
	}

	return 0.;

}

double I_1_analytic(double alpha, double delta)
{
	double A     = sqrt(alpha*alpha + delta*delta);
	double minus = sqrt(A - alpha);
	double plus  = sqrt(A + alpha);

	if(alpha == 0. || delta == 0. || minus == 0.)
	{
		return 0.;
	}

	double ans = 2. * ( (2. * alpha*alpha + (alpha*alpha - 1.)*delta*delta + pow(delta, 4.))*sin(A) 
                - (2. * alpha*alpha - delta*delta) * A * cos(A)) / pow(A, 5.);
    	return ans;

}

double I_1_of_2(double alpha, double delta)
{
	double A   = sqrt(alpha*alpha + delta*delta);
	double ans = -2. * delta*delta * (3. * A * cos(A) + (-3. + A*A) * sin(A)) / pow(A, 5.);
	return ans;
}

double I_2_analytic(double alpha, double delta)
{
	double A     = sqrt(alpha*alpha + delta*delta);
	double plus  = sqrt(A + alpha);
	double minus = sqrt(A - alpha);

	if(alpha == 0. || delta == 0. || minus == 0.)
        {
                return 0.;
        }

	double num   = 2. * alpha * delta * (3. * A * cos(A) + (-3. + A*A) * sin(A));
	double denom = pow(A, 5.);
	double ans   = num / denom;
	return ans;
}

double I_3_analytic(double alpha, double delta)
{
	if(alpha == 0. || delta == 0.)
	{
		return 0.;
	}	

	double A        = sqrt(alpha*alpha + delta*delta);
	double term1    = 6. * alpha*alpha * cos(A) / pow(A, 4.);
	double term2    = -2. * cos(A) / (A*A);
	double term3    = 6. * delta*delta * sin(A) / pow(A, 5.);
	double term4    = -4. * sin(A) / pow(A, 3.);
	double term5    = 2. * alpha*alpha * sin(A) / pow(A, 3.);
	double ans      = term1 + term2 + term3 + term4 + term5;
	return ans;
}

double I_3_limit(double alpha, double delta)
{
	if(alpha == 0.)
	{
		return 0.;
	}

        double A        = fabs(alpha);
        double term5    = 2. * alpha*alpha * sin(A) / pow(A, 3.);
        double ans      = term5;
        return ans;
}

double chi_11_integrand(double tau_prime, void * parameters)
{
	struct params * params = (struct params*) parameters;

	double prefactor  = 1.; //should be 1j
	double beta       = sqrt(1. - pow(params->gamma, -2.));
	double alpha      = beta * cos(params->theta) * tau_prime * params->gamma;
	double delta      = 2. * params->omega/(params->epsilon * params->omega_c) 
			   * sin(params->theta) * params->gamma * beta 
			   * sin((params->epsilon * params->omega_c / params->omega) * tau_prime / (2.));

	double gamma_term = beta*beta * params->gamma * Df(params);
//	double tau_term   = exp(1j * tau_prime * gamma) * sin((epsilon * omega_c / omega) * tau_prime);
//	double tau_term   = -sin(tau_prime * params->gamma) 
//			    * sin((epsilon * params->omega_c / params->omega) * tau_prime);
	double tauxi_term = 0.5 * (cos((params->epsilon * params->omega_c / params->omega) * tau_prime) * I_1_analytic(alpha, delta)
				   - I_1_of_2(alpha, delta));
	double ans        = prefactor * gamma_term * tauxi_term * params->gamma*params->gamma * beta;
	
	return ans;
}

double chi_12_integrand(double tau_prime, void * parameters)
{
	struct params * params = (struct params*) parameters;

	double prefactor  = 1.; //should be 1j
	double beta       = sqrt(1. - pow(params->gamma, -2.));
	double alpha      = beta * cos(params->theta) * tau_prime * params->gamma;
	double delta      = 2. * params->omega/(params->epsilon * params->omega_c) 
			   * sin(params->theta) * params->gamma * beta 
			   * sin((params->epsilon * params->omega_c / params->omega) * tau_prime / (2.));

	double gamma_term = beta*beta * params->gamma * Df(params);
//	double tau_term   = exp(1j * tau_prime * gamma) * sin((epsilon * omega_c / omega) * tau_prime);
//	double tau_term   = -sin(tau_prime * params->gamma) 
//			    * sin((epsilon * params->omega_c / params->omega) * tau_prime);
	double tau_term   = sin((params->epsilon * params->omega_c / params->omega) * tau_prime);
	double xi_term    = -0.5 * I_1_analytic(alpha, delta);
	double ans        = prefactor * gamma_term * xi_term * tau_term * params->gamma*params->gamma * beta;
	
	return ans;
}

double chi_13_integrand(double tau_prime, void * parameters)
{
	struct params * params = (struct params*) parameters;

	double prefactor  = 1.; //should be 1j
	double beta       = sqrt(1. - pow(params->gamma, -2.));
	double alpha      = beta * cos(params->theta) * tau_prime * params->gamma;
	double delta      = 2. * params->omega/(params->epsilon * params->omega_c) 
			   * sin(params->theta) * params->gamma * beta 
			   * sin((params->epsilon * params->omega_c / params->omega) * tau_prime / (2.));

	double gamma_term = beta*beta * params->gamma * Df(params);//* exp(-params->gamma/params->theta_e);
				//explicit imag part
//	double tau_term   = cos(tau_prime * params->gamma) * cos((params->epsilon * params->omega_c / params->omega) * tau_prime/2.);

//	double tau_term   = exp(1j * tau_prime * gamma) * sin((epsilon * omega_c / omega) * tau_prime);
//	double tau_term   = sin(tau_prime * params->gamma) 
//			    * sin((params->epsilon * params->omega_c / params->omega) * tau_prime / 2.);

	/*K_13 is K_32 except -S_2 -> C_2*/
	double tau_term   = cos((params->epsilon * params->omega_c / params->omega) * tau_prime / 2.);
	double xi_term    = -I_2_analytic(alpha, delta); //this term has a 1j, and I_2 has a 1j = -1 factor 
	double ans        = prefactor * gamma_term * xi_term * tau_term * params->gamma*params->gamma * beta;

	return ans;
}

/*note: for imaginary part of chi_22 splitting the integrand up is faster (40sec vs 170sec)
	but for real part the combined integrand is faster (3sec vs 10sec)*/
double chi_22_integrand_real(double tau_prime, void * parameters)
{
	struct params * params = (struct params*) parameters;

	double prefactor  = 1.; //should be 1j
	double beta       = sqrt(1. - pow(params->gamma, -2.));
	double alpha      = beta * cos(params->theta) * tau_prime * params->gamma;
	double delta      = 2. * params->omega/(params->epsilon * params->omega_c) 
			   * sin(params->theta) * params->gamma * beta 
			   * sin((params->epsilon * params->omega_c / params->omega) * tau_prime / (2.));

	double gamma_term = beta*beta * params->gamma * Df(params);
//	double tau_term   = exp(1j * tau_prime * gamma) * sin((epsilon * omega_c / omega) * tau_prime);
//	double tau_term   = -sin(tau_prime * params->gamma) 
//			    * sin((epsilon * params->omega_c / params->omega) * tau_prime);
	double tauxi_term = 0.5 * (cos((params->epsilon * params->omega_c / params->omega) * tau_prime) * I_1_analytic(alpha, delta)
				   + I_1_of_2(alpha, delta));

	double ans        = prefactor * gamma_term * tauxi_term * params->gamma*params->gamma * beta;
	
	return ans;
}

double chi_22_integrand_p1(double tau_prime, void * parameters)
{
        struct params * params = (struct params*) parameters;

        double prefactor  = 1.; //should be 1j
        double beta       = sqrt(1. - pow(params->gamma, -2.));
        double alpha      = beta * cos(params->theta) * tau_prime * params->gamma;
        double delta      = 2. * params->omega/(params->epsilon * params->omega_c)
                           * sin(params->theta) * params->gamma * beta 
                           * sin((params->epsilon * params->omega_c / params->omega) * tau_prime / (2.));

        double gamma_term = beta*beta * params->gamma * MJ(params);
//      double tau_term   = exp(1j * tau_prime * gamma) * sin((epsilon * omega_c / omega) * tau_prime);
//      double tau_term   = -sin(tau_prime * params->gamma) 
//                          * sin((epsilon * params->omega_c / params->omega) * tau_prime);
        double tauxi_term = 0.5 * (cos((params->epsilon * params->omega_c / params->omega) * tau_prime) * I_1_analytic(alpha, delta)
                                   + 0.);

        double ans        = prefactor * gamma_term * tauxi_term * params->gamma*params->gamma * beta;

        return ans;
}

double chi_22_integrand_p2(double tau_prime, void * parameters)
{
        struct params * params = (struct params*) parameters;

        double prefactor  = 1.; //should be 1j
        double beta       = sqrt(1. - pow(params->gamma, -2.));
        double alpha      = beta * cos(params->theta) * tau_prime * params->gamma;
        double delta      = 2. * params->omega/(params->epsilon * params->omega_c)
                           * sin(params->theta) * params->gamma * beta 
                           * sin((params->epsilon * params->omega_c / params->omega) * tau_prime / (2.));

        double gamma_term = beta*beta * params->gamma * MJ(params);
//      double tau_term   = exp(1j * tau_prime * gamma) * sin((epsilon * omega_c / omega) * tau_prime);
//      double tau_term   = -sin(tau_prime * params->gamma) 
//                          * sin((epsilon * params->omega_c / params->omega) * tau_prime);
        double tauxi_term = 0.5 * (0.
                                   + I_1_of_2(alpha, delta));

        double ans        = prefactor * gamma_term * tauxi_term * params->gamma*params->gamma * beta;

        return ans;
}

double chi_32_integrand(double tau_prime, void * parameters)
{
	struct params * params = (struct params*) parameters;

	double prefactor  = 1.; // should be 1j 
	double beta       = sqrt(1. - pow(params->gamma, -2.));
	double alpha      = beta * cos(params->theta) * tau_prime * params->gamma;
	double delta      = 2. * params->omega/(params->epsilon * params->omega_c) 
			   * sin(params->theta) * params->gamma * beta 
			   * sin((params->epsilon * params->omega_c / params->omega) * tau_prime / (2.));

	double gamma_term = beta*beta * params->gamma * Df(params);//* exp(-params->gamma/params->theta_e);
//	double tau_term   = exp(1j * tau_prime * gamma) * sin((epsilon * omega_c / omega) * tau_prime);
//	double tau_term   = sin(tau_prime * params->gamma) 
//			    * sin((params->epsilon * params->omega_c / params->omega) * tau_prime / 2.);
	double tau_term   = sin((params->epsilon * params->omega_c / params->omega) * tau_prime / 2.);
	double xi_term    = I_2_analytic(alpha, delta); //should be times 1j * -1j = +1
	double ans        = prefactor * gamma_term * xi_term * tau_term * params->gamma*params->gamma * beta;

	return ans;
}

double chi_33_integrand(double tau_prime, void * parameters)
{
	struct params * params = (struct params*) parameters;

	double prefactor  = 1.; //should be 1j
	double beta       = sqrt(1. - pow(params->gamma, -2.));
	double alpha      = beta * cos(params->theta) * tau_prime * params->gamma;
	double delta      = 2. * params->omega/(params->epsilon * params->omega_c) 
			   * sin(params->theta) * params->gamma * beta 
			   * sin((params->epsilon * params->omega_c / params->omega) * tau_prime / (2.));

	double gamma_term = beta*beta * params->gamma * Df(params);
//	double tau_term   = exp(1j * tau_prime * gamma) * sin((epsilon * omega_c / omega) * tau_prime);
//	double tau_term   = -sin(tau_prime * params->gamma) 
//			    * sin((epsilon * params->omega_c / params->omega) * tau_prime);
	double tau_term   = 1.;
//	double xi_term    = I_3_analytic(alpha, delta) - I_3_limit(alpha, delta);
	double xi_term;

	/*subtracting off this term on the imaginary part speeds convergence
	  by a factor of 2-3.  The term integrates to zero, so we get this
	  speed boost basically for free.  The real part of this term is
	  nonzero, however, and actually slows convergence by a factor of 10.*/
	if(params->real == 1)
	{
		xi_term = I_3_analytic(alpha, delta);
		xi_term *= -sin(params->gamma * tau_prime);
	}
	else
	{
		xi_term = I_3_analytic(alpha, delta) - I_3_limit(alpha, delta);
		xi_term *= cos(params->gamma * tau_prime);
	}

	double ans        = prefactor * gamma_term * xi_term * tau_term * params->gamma*params->gamma * beta;
	
	return ans;
}
