#include <stdio.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_errno.h>
#include "susceptibility_tensor.h"

double tau_integrator(double gamma, void * parameters)
{
	struct params * params = (struct params*) parameters;

	if(gamma == 1.)
	{
		return 0.;
	}

        double ans_tot  = 0.;
	double ans_step = 0.;
	double error    = 0.;
	double step;
	double start    = 0.;
        size_t n        = 50;
        size_t limit    = 5000;
        double epsabs   = 0.;
        double epsrel   = 1e-8;
	double sign_correction;
	enum gsl_integration_qawo_enum gsl_weight;

	gsl_integration_workspace * w = gsl_integration_workspace_alloc (5000);

	if(params->integrand == &chi_33_integrand)
	{
		step = 2. * M_PI / gamma;
		sign_correction = 1.;
	}

	else
	{
		step     = M_PI/gamma;
	}

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


	//need to update value of gamma
	params-> gamma = gamma;

	gsl_set_error_handler_off();
	gsl_function F;
	F.function = params->integrand;
	F.params   = params;

	int i            = 0;
	int max_counter  = 1000;
	double tolerance = 1e-6;
	int counts       = 0;

	/*TODO: explain this*/
	double ans_sign         = 0; 
	int sign_change_counter = 0;
	int max_sign_change_counter = 1200.;

	int i_max        = 1000;
	double small_tol = 1e-17;
	while(i == 0 || counts < max_counter)
	{

		if(params->integrand == &chi_33_integrand)
		{
			gsl_integration_qng(&F, i*step, (i+1)*step, epsabs, epsrel, &ans_step, &error, &limit);
		}
		else
		{
			gsl_integration_qawo(&F, i*step, epsabs, epsrel, limit, w, table, &ans_step, &error);
		}
		ans_tot += ans_step;
		i       += 1;


		if(fabs(ans_step / ans_tot) < tolerance)
		{
			counts += 1;
		}

		if(i == 1 || ans_sign != ans_tot/fabs(ans_tot))
		{
			ans_sign = ans_tot/fabs(ans_tot);
			sign_change_counter++;
		}
		if(sign_change_counter >= max_sign_change_counter)
		{
			return 0.;
		}


	}

	gsl_integration_qawo_table_free(table);
	gsl_integration_workspace_free(w);
	ans_tot = ans_tot * sign_correction;

	return ans_tot;
}

double gamma_integrator(struct params *p)
{
        double prefactor = 2. * M_PI * p->omega_p*p->omega_p / (p->omega * p->omega);

        gsl_function F;
        F.function = &tau_integrator;
        F.params   = p;
	
//	gsl_integration_workspace * w = gsl_integration_workspace_alloc(5000);


        double start  = 1.;
        double end    = 150.; //figure out better way to do this
        double ans    = 0.;
        double error  = 0.;
        size_t limit  = 5000;
        double epsabs = 0.;
        double epsrel = 1e-3;


//        gsl_set_error_handler_off();

        gsl_integration_qng(&F, start, end, epsabs, epsrel, &ans, &error, &limit);

//      gsl_integration_qagiu(&F, start, epsabs, epsrel, limit, w, &ans, &error);
//      gsl_integration_workspace_free(w);

        return prefactor * ans;
}


double chi_11(struct params * p)
{
	p->integrand = &chi_11_integrand;

	return gamma_integrator(p);
}

double chi_12(struct params * p)
{
        p->integrand = &chi_12_integrand;

        return gamma_integrator(p);
}

double chi_32(struct params * p)
{
        p->integrand = &chi_32_integrand;

        return gamma_integrator(p);
}

double chi_13(struct params * p)
{
        p->integrand = &chi_13_integrand;

        return gamma_integrator(p);
}

double chi_22(struct params * p)
{
	double ans = 0.;

	if(p->real == 1)
	{
		p->integrand = &chi_22_integrand_p1;
		ans = gamma_integrator(p);
		p->integrand = &chi_22_integrand_p2;
		ans += gamma_integrator(p);
	}
	else
	{
		p->integrand = &chi_22_integrand_real;
		ans = gamma_integrator(p);
	}

	return ans;
}

//double chi_33(struct params * p)
//{
//        p->integrand = &chi_33_integrand;
//
//        return gamma_integrator(p);
//}

double chi_21(struct params * p)
{
        return -chi_12(p);
}

double chi_23(struct params * p)
{
        return -chi_32(p);
}

double chi_31(struct params * p)
{
        return chi_13(p);
}
//double chi_12(struct params * p)
//{
//	double prefactor = 2. * M_PI * p->omega_p*p->omega_p / (p->omega * p->omega);
//	
//	gsl_function F;
//	p->integrand = &chi_12_integrand;
//        F.function = &tau_integrator;
//        F.params   = p;
////	gsl_integration_workspace * w = gsl_integration_workspace_alloc(5000);
//
//
////	double start  = start_search_12(p);
//        double start  = 1.;
//	double end    = 150.; //figure out better way to do this
//	double ans    = 0.;
//	double error  = 0.;
//	size_t limit  = 50;
//	double epsabs = 0.;
//	double epsrel = 1e-8;
//
//	
////	gsl_set_error_handler_off();
//
//	gsl_integration_qng(&F, start, end, epsabs, epsrel, &ans, &error, &limit);
//
////	gsl_integration_qagiu(&F, start, 0., 1e-8, limit, w, &ans, &error);
////	gsl_integration_workspace_free(w);
//
//	return prefactor * ans;
//}
//
//double chi_21(struct params * p)
//{
//	return -chi_12(p);
//}
//
//double chi_13(struct params * p)
//{
//	double prefactor = 2. * M_PI * p->omega_p*p->omega_p / (p->omega * p->omega);
//	
//	gsl_function F;
//        F.function = &tau_integrator;
//	p->integrand = &chi_13_integrand;
//        F.params   = p;
////	gsl_integration_workspace * w = gsl_integration_workspace_alloc(5000);
//
//        double start  = 1.;
////	double start  = start_search_13(p);
//	double end    = 150.;
//	double ans    = 0.;
//	double error  = 0.;
//	size_t limit  = 50;
//	double epsabs = 0.;
//        double epsrel = 1e-8;
//
//
//	gsl_set_error_handler_off();
//	gsl_integration_qng(&F, start, end, epsabs, epsrel, &ans, &error, &limit);
//
////	gsl_integration_qagiu(&F, start, 0., 1e-8, limit, w, &ans, &error);
////	gsl_integration_workspace_free(w);
//
//	return prefactor * ans;
////	return ans;
//}
//
//double chi_31(struct params * p)
//{
//	return chi_13(p);
//}

//double tau_integrator_22_p1(double gamma, void * parameters)
//{
//        struct params * params = (struct params*) parameters;
//
//        if(gamma == 1.)
//        {
//                return 0.;
//        }
//
//
//        double ans_tot  = 0.;
//        double ans_step = 0.;
//        double error    = 0.;
//        double step     = M_PI / gamma; //TODO: change or play with this parameter
//        double start    = 0.;
////        double end      = M_PI * params->omega / params->omega_c * 2. * params->resolution_factor;
//        size_t n        = 50;
//        size_t limit    = 5000;
//        double epsabs   = 0.;
//        double epsrel   = 1e-8;
//        enum gsl_integration_qawo_enum gsl_weight;
//        double sign_correction;
//        //need to update value of gamma
//        params-> gamma = gamma;
//
//        /*set up GSL QAWO integrator.  Do we need a new table w every call to tau_integrator_12?*/
//        /*we should also try QAWF; it fits the integrals we need, and may be faster than QAWO.  */
//
//        if(params->real == 1)
//        {
//                gsl_weight      = GSL_INTEG_SINE;
//                sign_correction = -1.;
//        }
//        else
//        {
//                gsl_weight      = GSL_INTEG_COSINE;
//                sign_correction = 1.;
//        }
//
//        gsl_integration_qawo_table * table =
//                                gsl_integration_qawo_table_alloc(gamma, step, gsl_weight, n);
//
//        gsl_integration_workspace * w = gsl_integration_workspace_alloc (5000);
//        gsl_set_error_handler_off();
//        gsl_function F;
//        F.function = &chi_22_integrand_p1;
//        F.params   = params;
//
//        int i            = 0;
//        int max_counter  = 2000;
//        double tolerance = 1e-6;
//        int counts       = 0;
//
//	/*TODO: explain this*/
//        double ans_sign         = 0;
//        int sign_change_counter = 0;
//        int max_sign_change_counter = 10000.;
//
//        int i_max        = 1000;
//        double small_tol = 1e-17;
//        while(i == 0 || counts < max_counter)
//        {
//                gsl_integration_qawo(&F, i*step, epsabs, epsrel, limit, w, table, &ans_step, &error);
//                ans_tot += ans_step;
//                i       += 1;
//
//                if(fabs(ans_step / ans_tot) < tolerance)
//                {
//                        counts += 1;
//                }
//
////                if(i >= i_max && fabs(ans_tot) < small_tol)
////                {
////                        counts = max_counter;
////                }
//
//		if(i == 1 || ans_sign != ans_tot/fabs(ans_tot))
//                {
//                        ans_sign = ans_tot/fabs(ans_tot);
//                        sign_change_counter++;
//                }
//                if(sign_change_counter >= max_sign_change_counter)
//                {
//                        return 0.;
//                }
//
//        }
//
//        gsl_integration_qawo_table_free(table);
//        gsl_integration_workspace_free(w);
//
//        return ans_tot * sign_correction;
//}
//
//double tau_integrator_22_p2(double gamma, void * parameters)
//{
//        struct params * params = (struct params*) parameters;
//
//        if(gamma == 1.)
//        {
//                return 0.;
//        }
//
//
//        double ans_tot  = 0.;
//        double ans_step = 0.;
//        double error    = 0.;
//        double step     = M_PI / gamma; //TODO: change or play with this parameter
//        double start    = 0.;
////        double end      = M_PI * params->omega / params->omega_c * 2. * params->resolution_factor;
//        size_t n        = 50;
//        size_t limit    = 5000;
//        double epsabs   = 0.;
//        double epsrel   = 1e-8;
//        enum gsl_integration_qawo_enum gsl_weight;
//        double sign_correction;
//        //need to update value of gamma
//        params-> gamma = gamma;
//
//        /*set up GSL QAWO integrator.  Do we need a new table w every call to tau_integrator_12?*/
//        /*we should also try QAWF; it fits the integrals we need, and may be faster than QAWO.  */
//
//        if(params->real == 1)
//        {
//                gsl_weight      = GSL_INTEG_SINE;
//                sign_correction = -1.;
//        }
//        else
//        {
//                gsl_weight      = GSL_INTEG_COSINE;
//                sign_correction = 1.;
//        }
//
//        gsl_integration_qawo_table * table =
//                                gsl_integration_qawo_table_alloc(gamma, step, gsl_weight, n);
//
//        gsl_integration_workspace * w = gsl_integration_workspace_alloc (5000);
//        gsl_set_error_handler_off();
//        gsl_function F;
//        F.function = &chi_22_integrand_p2;
//        F.params   = params;
//
//        int i            = 0;
//        int max_counter  = 2000;
//        double tolerance = 1e-6;
//        int counts       = 0;
//
//	/*TODO: explain this*/
//        double ans_sign         = 0;
//        int sign_change_counter = 0;
//        int max_sign_change_counter = 10000.;
//
//        int i_max        = 1000;
//        double small_tol = 1e-17;
//        while(i == 0 || counts < max_counter)
//        {
//                gsl_integration_qawo(&F, i*step, epsabs, epsrel, limit, w, table, &ans_step, &error);
//                ans_tot += ans_step;
//                i       += 1;
//
//                if(fabs(ans_step / ans_tot) < tolerance)
//                {
//                        counts += 1;
//                }
//
////                if(i >= i_max && fabs(ans_tot) < small_tol)
////                {
////                        counts = max_counter;
////                }
//        
//		if(i == 1 || ans_sign != ans_tot/fabs(ans_tot))
//                {
//                        ans_sign = ans_tot/fabs(ans_tot);
//                        sign_change_counter++;
//                }
//                if(sign_change_counter >= max_sign_change_counter)
//                {
//                        return 0.;
//                }
//
//	}
//
//        gsl_integration_qawo_table_free(table);
//        gsl_integration_workspace_free(w);
//
//        return ans_tot * sign_correction;
//}
//
//double tau_integrator_22_real(double gamma, void * parameters)
//{
//	struct params * params = (struct params*) parameters;
//
//	if(gamma == 1.)
//	{
//		return 0.;
//	}
//
//
//        double ans_tot  = 0.;
//	double ans_step = 0.;
//	double error    = 0.;
//        double step     = M_PI / gamma; //TODO: change or play with this parameter
//        double start    = 0.;
////        double end      = M_PI * params->omega / params->omega_c * 2. * params->resolution_factor;
//	size_t n        = 50;
//	size_t limit    = 5000;
//	double epsabs   = 0.;
//	double epsrel   = 1e-8;
//	enum gsl_integration_qawo_enum gsl_weight;
//	double sign_correction;
//	//need to update value of gamma
//	params-> gamma = gamma;
//
//	/*set up GSL QAWO integrator.  Do we need a new table w every call to tau_integrator_12?*/
//	/*we should also try QAWF; it fits the integrals we need, and may be faster than QAWO.  */
//
//	if(params->real == 1)
//	{
//		gsl_weight      = GSL_INTEG_SINE;
//		sign_correction = -1.;
//	}
//	else
//	{
//		gsl_weight      = GSL_INTEG_COSINE;
//		sign_correction = 1.;
//	}
//
//	gsl_integration_qawo_table * table = 
//				gsl_integration_qawo_table_alloc(gamma, step, gsl_weight, n);
//	
//	gsl_integration_workspace * w = gsl_integration_workspace_alloc (5000);
//	gsl_set_error_handler_off();
//	gsl_function F;
//	F.function = &chi_22_integrand_real;
//	F.params   = params;
//
//	int i            = 0;
//	int max_counter  = 1000;
//	double tolerance = 1e-6;
//	int counts       = 0;
//
//	/*TODO: explain this*/
//        double ans_sign         = 0;
//        int sign_change_counter = 0;
//        int max_sign_change_counter = 10000.;
//
//	int i_max        = 1000;
//	double small_tol = 1e-20;
//	while(i == 0 || counts < max_counter)
//	{
//		gsl_integration_qawo(&F, i*step, epsabs, epsrel, limit, w, table, &ans_step, &error);
//		ans_tot += ans_step;
//		i       += 1;
//
//		if(fabs(ans_step / ans_tot) < tolerance)
//		{
//			counts += 1;
//		}
//
//		if(i >= i_max && fabs(ans_tot) < small_tol)
//		{
//			counts = max_counter;
//		}
//
////		if(i == 1 || ans_sign != ans_tot/fabs(ans_tot))
////                {
////                        ans_sign = ans_tot/fabs(ans_tot);
////                        sign_change_counter++;
////                }
////                if(sign_change_counter >= max_sign_change_counter)
////                {
////                        return 0.;
////                }
//	}
//
//	gsl_integration_qawo_table_free(table);
//	gsl_integration_workspace_free(w);
//
//	return ans_tot * sign_correction;
//}
//
//double tau_integrator_22(double gamma, void * parameters)
//{
//	struct params * params = (struct params*) parameters;
//
//	double ans = 0.;
//
//	if(params->real == 0)
//	{
//		params->integrand = &chi_22_integrand_p1;
//		ans = tau_integrator(gamma, &params);
//
//		params->integrand = &chi_22_integrand_p2;
//		ans += tau_integrator(gamma, &params);
//	
////		ans = tau_integrator_22_p1(gamma, parameters)
////                     +tau_integrator_22_p2(gamma, parameters);
//	}
//	else
//	{
//		params->integrand = &chi_22_integrand_real;
//		ans = tau_integrator(gamma, &params);
////		ans = tau_integrator_22_real(gamma, parameters);
//	}
//	
//	return ans;
//}
//
//double chi_22(struct params * p)
//{
//	double prefactor = 2. * M_PI * p->omega_p*p->omega_p / (p->omega * p->omega);
//	
//	gsl_function F;
//        F.function = &tau_integrator_22;
//        F.params   = p;
////	gsl_integration_workspace * w = gsl_integration_workspace_alloc(5000);
//
//
//	double start  = 1.;
////	double start  = start_search_22(p);
//	double end    = 150.; //figure out better way to do this
//	double ans    = 0.;
//	double error  = 0.;
//	size_t limit  = 50;
//	double epsabs = 0.;
//	double epsrel = 1e-8;
//
//	
//	gsl_set_error_handler_off();
//
//	gsl_integration_qng(&F, start, end, epsabs, epsrel, &ans, &error, &limit);
//
////	gsl_integration_qagiu(&F, start, 0., 1e-8, limit, w, &ans, &error);
////	gsl_integration_workspace_free(w);
//
//	return prefactor * ans;
//}

//double chi_32(struct params * p)
//{
//	double prefactor = 2. * M_PI * p->omega_p*p->omega_p / (p->omega * p->omega);
//	
//	gsl_function F;
//	p->integrand = &chi_32_integrand;
//        F.function = &tau_integrator;
//        F.params   = p;
////	gsl_integration_workspace * w = gsl_integration_workspace_alloc(5000);
//
////	double start  = start_search_32(p);
//        double start  = 1.;
//	double end    = 150.;
//	double ans    = 0.;
//	double error  = 0.;
//	size_t limit  = 50;
//	double epsabs = 0.;
//        double epsrel = 1e-8;
//
////	gsl_set_error_handler_off();
//
//	gsl_integration_qng(&F, start, end, epsabs, epsrel, &ans, &error, &limit);
//
////	gsl_integration_qagiu(&F, start, 0., 1e-8, limit, w, &ans, &error);
////	gsl_integration_workspace_free(w);
//
//	return prefactor * ans;
//}
//
//double chi_23(struct params * p)
//{
//	return -chi_32(p);
//}
//
double chi_33(struct params * p)
{
	double prefactor = 2. * M_PI * p->omega_p*p->omega_p / (p->omega * p->omega);
	
	gsl_function F;
        F.function = &tau_integrator;
	p->integrand = &chi_33_integrand;
        F.params   = p;
//	gsl_integration_workspace * w = gsl_integration_workspace_alloc(5000);


	double start  = 1.; //start_search_12(p);
	double end    = 150.; //figure out better way to do this
	double ans    = 0.;
	double error  = 0.;
	size_t limit  = 50;
	double epsabs = 0.;
	double epsrel = 1e-8;

	
	gsl_set_error_handler_off();

	gsl_integration_qng(&F, start, end, epsabs, epsrel, &ans, &error, &limit);

//	gsl_integration_qagiu(&F, start, 0., 1e-8, limit, w, &ans, &error);
//	gsl_integration_workspace_free(w);

	return prefactor * ans;
}
