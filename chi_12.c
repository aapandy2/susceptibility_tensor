#include <stdio.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_errno.h>
#include "susceptibility_tensor.h"

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

double tau_integrator_12(double gamma, void * parameters)
{
	struct params * params = (struct params*) parameters;

	if(gamma == 1.)
	{
		return 0.;
	}

        double ans_tot  = 0.;
	double ans_step = 0.;
	double error    = 0.;
        double step     = M_PI/gamma; //TODO: change or play with this parameter
        double start    = 0.;
//        double end      = M_PI * params->omega / params->omega_c * 2. * params->resolution_factor;
	size_t n        = 50;
	size_t limit    = 5000;
	double epsabs   = 0.;
	double epsrel   = 1e-8;
	enum gsl_integration_qawo_enum gsl_weight;
	double sign_correction;
	//need to update value of gamma
	params-> gamma = gamma;

	/*set up GSL QAWO integrator.  Do we need a new table w every call to tau_integrator_12?*/
	/*we should also try QAWF; it fits the integrals we need, and may be faster than QAWO.  */

	if(params->real == 1)
	{
		gsl_weight      = GSL_INTEG_SINE;
		sign_correction = -1.;
	}
	else
	{
		gsl_weight      = GSL_INTEG_COSINE;
		sign_correction = 1.;
	}

	gsl_integration_qawo_table * table = 
				gsl_integration_qawo_table_alloc(gamma, step, gsl_weight, n);
	
	gsl_integration_workspace * w = gsl_integration_workspace_alloc (5000);
	gsl_set_error_handler_off();
	gsl_function F;
	F.function = &chi_12_integrand;
	F.params   = params;

	int i            = 0;
	int max_counter  = 500;
	double tolerance = 1e-5;
	int counts       = 0;

	int i_max        = 1000;
	double small_tol = 1e-20;
	while(i == 0 || counts < max_counter)
	{
		gsl_integration_qawo(&F, i*step, epsabs, epsrel, limit, w, table, &ans_step, &error);
		ans_tot += ans_step;
		i       += 1;

		if(fabs(ans_step / ans_tot) < tolerance)
		{
			counts += 1;
		}

		if(i >= i_max && fabs(ans_tot) < small_tol)
		{
			counts = max_counter;
		}
	}


	gsl_integration_qawo_table_free(table);
	gsl_integration_workspace_free(w);

	return ans_tot * sign_correction;
}

double chi_12(struct params * p)
{
	double prefactor = 2. * M_PI * p->omega_p*p->omega_p / (p->omega * p->omega);
	
	gsl_function F;
        F.function = &tau_integrator_12;
        F.params   = p;
//	gsl_integration_workspace * w = gsl_integration_workspace_alloc(5000);


//	double start  = start_search_12(p);
        double start  = 1.;
	double end    = 10.; //figure out better way to do this
	double ans    = 0.;
	double error  = 0.;
	size_t limit  = 50;
	double epsabs = 0.;
	double epsrel = 1e-8;

	
//	gsl_set_error_handler_off();

	gsl_integration_qng(&F, start, end, epsabs, epsrel, &ans, &error, &limit);

//	gsl_integration_qagiu(&F, start, 0., 1e-8, limit, w, &ans, &error);
//	gsl_integration_workspace_free(w);

	return prefactor * ans;
}

double chi_21(struct params * p)
{
	return -chi_12(p);
}
