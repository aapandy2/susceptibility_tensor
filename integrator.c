#include <stdio.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_errno.h>
#include "susceptibility_tensor.h"
# include <omp.h>

double gamma_integrand_gsl(double gamma, void * parameters)
{
  struct params * params = (struct params*) parameters;

  double omega = params->omega;

  return gamma_integrand(gamma, omega, parameters);
}

/*gamma_integrand: integrator for the first integral (the tau integral) for the
 *                components of the susceptibility tensor.  The chi_33 integral
 *                is faster when we use a fixed-order Gaussian quadrature
 *                method GSL QNG, rather than the GSL QAWO adaptive integrator
 *                used for the other components.
 *
 *@params: double gamma (the second integration is carried out over this
 *         variable), void * parameters (a pointer to the struct of parameters)
 *
 *@returns: numerically evaluated tau integral at a given gamma of a given
 *          component of the susceptibility tensor
 */
double gamma_integrand(double gamma, double omega, void * parameters)
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
  
  if(params->tau_integrand == &chi_33_integrand)
  {
    step = 2. * M_PI / gamma;
    sign_correction = 1.;
  }
  else
  {
    step     = M_PI/gamma;
  }
  
  gsl_integration_qawo_table * table =
                   gsl_integration_qawo_table_alloc(gamma, step, gsl_weight, n);
  
  
//need to update value of gamma
//	params-> gamma = gamma;
  struct params tau_params = *params;
  tau_params.gamma = gamma;
  tau_params.omega = omega;
  
  gsl_set_error_handler_off();
  gsl_function F;
  F.function = params->tau_integrand;
  F.params   = &tau_params;
  
  int i            = 0;
  int max_counter  = 2500;
  double tolerance = 1e-7;
  int counts       = 0;
  
  /*TODO: explain this */
  int small_counter = 0;
  double small_tol  = 1e-7;  //TODO: want to set this adaptively
  int max_small_counter = 1000;
  
  while(i == 0 || counts < max_counter)
  {
  
    if(params->tau_integrand == &chi_33_integrand)
    {
      gsl_integration_qng(&F, i*step, (i+1)*step, epsabs, epsrel, &ans_step, 
                          &error, &limit);
    }
    else
    {
      gsl_integration_qawo(&F, i*step, epsabs, epsrel, limit, w, table, 
                           &ans_step, &error);
    }
    
    ans_tot += ans_step;
    i       += 1;
    
    if(fabs(ans_step / ans_tot) < tolerance)
    {
      counts += 1;
    }
    
    if(fabs(ans_tot) < small_tol)
    {
      small_counter++;
    }
    if(small_counter >= max_small_counter)
    {
      return 0.;
    }
  }
  
  gsl_integration_qawo_table_free(table);
  gsl_integration_workspace_free(w);
  ans_tot = ans_tot * sign_correction;
  
  return ans_tot;
}

double midpoint_rule(struct params *params, double start, 
                     double end, int samples)
{
  double a = start;
  double b = end;
  int i;
  int n = samples;
  double total;
  double x;
  total = 0.0;

  int tid;

# pragma omp parallel for private(i , x) reduction ( + : total )

  for(i = 0; i < n; i++ )
  {
    x = (b - a)/n * i + a;
    total = total + gamma_integrand(x, params->omega, params);
  }

  total = (b - a) / n * total;

  return total;
}

//TODO: can we define NUM_QUAD without a hash define?
#define NUM_QUAD 21
double gauss_legendre(struct params * params, double start, double end)
{
  double quadPts[NUM_QUAD] = \
        {-9.93752171e-01,  -9.67226839e-01,  -9.20099334e-01,
         -8.53363365e-01,  -7.68439963e-01,  -6.67138804e-01,
         -5.51618836e-01,  -4.24342120e-01,  -2.88021317e-01,
         -1.45561854e-01,   1.98918497e-16,   1.45561854e-01,
          2.88021317e-01,   4.24342120e-01,   5.51618836e-01,
          6.67138804e-01,   7.68439963e-01,   8.53363365e-01,
          9.20099334e-01,   9.67226839e-01,   9.93752171e-01};

  double weights[NUM_QUAD] = \
        {0.01601723,  0.03695379,  0.05713443,  0.07610011,  0.09344442,
         0.1087973 ,  0.12183142,  0.13226894,  0.13988739,  0.1445244 ,
         0.14608113,  0.1445244 ,  0.13988739,  0.13226894,  0.12183142,
         0.1087973 ,  0.09344442,  0.07610011,  0.05713443,  0.03695379,
         0.01601723};

  /*first we change integration variables to x = 1/(gamma^2 - 1), where the
    upper integration bound b = 1 and the lower bound a = 0.  We then apply
    the transformation: 
    \int_a^b f(x) dx = (b-a)/2 \int_{-1}^1 f((b-a)x/2 + (a+b)/2) dx
                     =     1/2 \int_{-1}^1 f(     x/2 +     1/2) dx */
  double x      = 0.;
  int i         = 0;
  double weight = 0.;
  double sum    = 0.;
  int n = NUM_QUAD;

  # pragma omp parallel for private(i , x, weight) reduction ( + : sum )

  for(i = 0; i < n; i++)
  {
    x        = quadPts[i];
    weight   = weights[i];

    sum = sum + (end - start)/2. * gamma_integrand((end - start)/2. * x + (end + start)/2., params->omega, params)
                      * weight;
  }

  return sum;
}

double simpsons_rule(struct params * params, double start, double end, 
                     int samples)
{
  int n,i;
  float s1=0,s2=0,sum,a,b,h;
  double upper_limit = end;
  double lower_limit = start;
  int num_points  = samples; //TODO: enforce that this value is even
  b = end;
  a = start;
  n = num_points;
  double total = 0.;
  
  int tid;
  
  h=(upper_limit-lower_limit)/num_points;
  
  if(num_points%2 != 0)
  {
    printf("the rule is not applicable");
  }
  
  # pragma omp parallel for private(i, s1, s2) reduction ( + : total )
  
  for(i=1;i<=num_points-1;i++)
  {
    if(i%2==0)
    {
      total = total + 2. * gamma_integrand(lower_limit + i*h, params->omega, params);
    }
    else
    {
      total = total + 4. * gamma_integrand(lower_limit + i*h, params->omega, params);
    }
  }

  sum = h/3 * (gamma_integrand(lower_limit, params->omega, params) 
        + gamma_integrand(upper_limit, params->omega, params) + total);

  return sum;
}

/*end_approx: approximate ending point for the gamma integral, beyond which
 *            contributions to the integral are assumed to be negligible.
 *
 *@params: pointer to struct of parameters *params
 *
 *@returns: approximate end location for gamma integral
 */
double end_approx(struct params *params)
{
  double end;
  
  double MJ_max        = 0.5 * (3. * params->theta_e 
                         + sqrt(4. + 9. * params->theta_e * params->theta_e));
  double PL_max_real   = sqrt((1. + params->pl_p)/params->pl_p);
  double PL_max_moving = 50./sqrt(params->omega/params->omega_c) 
                         + 9. * pow(params->omega/params->omega_c, 1./3.);
  double kappa_max     = (-3. + 3. * params->kappa_width * params->kappa 
  		          + sqrt(1. - 4. * params->kappa 
  		  	         - 18. * params->kappa_width * params->kappa 
                                 + 4. * pow(params->kappa, 2.) 
                                 + 9. * pow(params->kappa_width
                                            * params->kappa, 2.))) 
  		         / (2. * (params->kappa - 2.));
  
  if(params->dist == 0)
  {
    end = 7. * MJ_max;
  }
  else if(params->dist == 1)
  {
    end = PL_max_moving;
  }
  else if(params->dist == 2)
  {
    end = 7. * kappa_max;
  }
  else
  {
    printf("\ndistribution or real/imag is set incorrectly");
    return 1.;
  }
  
  return end;
}

/*gsl_integrator: wrapper for a GSL integrator, to be used for the gamma
 *                integral.  This function can be very expensive, because
 *                it makes a large number of function calls to a region
 *                of parameter space that has a slowly convergent tau integral.
 *
 *@params: pointer to struct of parameters *params, start (starting point for 
 *         gamma integral), end (ending point for gamma integral)
 *
 *@returns: gamma integral for a given component of chi_ij, evaluated using
 *          a GSL adaptive integrator
 */
double gsl_integrator(struct params *params, double start, double end)
{
  gsl_function F;
  F.function = gamma_integrand_gsl;
  F.params   = params;
  
  gsl_integration_workspace * w = gsl_integration_workspace_alloc(5000);
  
  double ans    = 0.;
  double error  = 0.;
  size_t limit  = 5000;
  double epsabs = 0.;
  double epsrel = 1e-5;
  int gsl_key   = 1;
  
  printf("\n%e\n", end);
  
  
  gsl_integration_qags(&F, start, end, epsabs, epsrel, limit, w, &ans, &error);

//gsl_integration_qng(&F, start, end, epsabs, epsrel, &ans, &error, &limit);
//gsl_integration_qag(&F, start, end, epsabs, epsrel, limit, gsl_key, w, &ans, &error);
//gsl_integration_qagiu(&F, start, epsabs, epsrel, limit, w, &ans, &error);
	
  gsl_integration_workspace_free(w);
  
  return ans;
}

/*gamma_integrator: function that evaluates the gamma integral, using one of
 *                  the above integration algorithms.
 *
 *@params: pointer to struct of parameters *params
 *
 *@returns: gamma integral for a given component of chi_ij
 */
double gamma_integrator(struct params *params)
{
  double prefactor = 2. * M_PI * params->omega_p*params->omega_p 
                     / (params->omega * params->omega);
  
  double start  = 1.;
  double end = end_approx(params);
  double ans_tot;
  
  /*for power-law, there is a singularity at low frequency at gamma = 1,
    which cannot be resolved by a trapezoidal integrator.  We are forced
    to use the (much slower) GSL integrator QAGS in these cases. */
  if(params->dist == 1 && params->omega/params->omega_c < 5.)
  {
    ans_tot = gsl_integrator(params, start, end);
  }
  else
  {
//    ans_tot = simpsons_rule(params, start, end, 150);
//		ans_tot = midpoint_rule(params, start, end, 150);
    ans_tot = gauss_legendre(params, start, end);
  }

//	double ans_tot = trapezoidal(params, start, end, 100);
//	double ans_tot = trapezoidal_adaptive(params, start, 1.);
//	double ans_tot = gsl_integrator(params, start, end);

  return prefactor * ans_tot;
}

/*chi_11: evaluates the component chi_11 of the susceptibility tensor.
 *
 *@params: pointer to struct of parameters *params
 *
 *@returns: the component chi_11, evaluated using the parameters in struct p
 */
double chi_11(struct params * params)
{
  if(params->use_spline == 1)
  {
  //  return chi_11_spline(params);
  }

  params->tau_integrand = &chi_11_integrand;
  params->gamma_integrand = &gamma_integrand;
  
  return gamma_integrator(params);
}

/*chi_12: evaluates the component chi_12 of the susceptibility tensor.
 *
 *@params: pointer to struct of parameters *params
 *
 *@returns: the component chi_12, evaluated using the parameters in struct p
 */
double chi_12(struct params * params)
{
  if(params->use_spline == 1)
  {
  //  return chi_12_spline(params);
  }

  params->tau_integrand = &chi_12_integrand;
  params->gamma_integrand = &gamma_integrand;
  
  return gamma_integrator(params);
}

/*chi_32: evaluates the component chi_32 of the susceptibility tensor.
 *
 *@params: pointer to struct of parameters *params
 *
 *@returns: the component chi_32, evaluated using the parameters in struct p
 */
double chi_32(struct params * params)
{
  if(params->use_spline == 1)
  {
  //  return chi_32_spline(params);
  }  

  params->tau_integrand = &chi_32_integrand;
  params->gamma_integrand = &gamma_integrand;
  
  return gamma_integrator(params);
}

/*chi_13: evaluates the component chi_13 of the susceptibility tensor.
 *
 *@params: pointer to struct of parameters *params
 *
 *@returns: the component chi_13, evaluated using the parameters in struct p
 */
double chi_13(struct params * params)
{
  if(params->use_spline == 1)
  {
  //  return chi_13_spline(params);
  }

  params->tau_integrand = &chi_13_integrand;
  params->gamma_integrand = &gamma_integrand;
  
  return gamma_integrator(params);
}

/*chi_22: evaluates the component chi_22 of the susceptibility tensor.
 *
 *@params: pointer to struct of parameters *params
 *
 *@returns: the component chi_22, evaluated using the parameters in struct p
 */
double chi_22(struct params * params)
{
  if(params->use_spline == 1)
  {
  //  return chi_22_spline(params);
  }

  double ans = 0.;
  params->gamma_integrand = &gamma_integrand;
  
  if(params->real == 0)
  {
    params->tau_integrand = &chi_22_integrand_p1;
    ans = gamma_integrator(params);
    params->tau_integrand = &chi_22_integrand_p2;
    ans += gamma_integrator(params);
  }
  else
  {
    params->tau_integrand = &chi_22_integrand_real;
    ans = gamma_integrator(params);
  }
  
  return ans;
}

/*chi_33: evaluates the component chi_33 of the susceptibility tensor.
 *
 *@params: pointer to struct of parameters *params
 *
 *@returns: the component chi_33, evaluated using the parameters in struct p
 */
double chi_33(struct params * params)
{
  if(params->use_spline == 1)
  {
  //  return chi_33_spline(params);
  }

  params->tau_integrand = &chi_33_integrand;
  params->gamma_integrand = &gamma_integrand;
  
  return gamma_integrator(params);
}

/*chi_21: evaluates the component chi_21 of the susceptibility tensor.
 *
 *@params: pointer to struct of parameters *params
 *
 *@returns: the component chi_21, evaluated using the parameters in struct p
 */
double chi_21(struct params * params)
{
  if(params->use_spline == 1)
  {
  //  return -chi_12_spline(params);
  }

  return -chi_12(params);
}

/*chi_23: evaluates the component chi_23 of the susceptibility tensor.
 *
 *@params: pointer to struct of parameters *params
 *
 *@returns: the component chi_23, evaluated using the parameters in struct p
 */
double chi_23(struct params * params)
{
  if(params->use_spline == 1)
  {
  //  return -chi_32_spline(params);
  }

  return -chi_32(params);
}

/*chi_31: evaluates the component chi_31 of the susceptibility tensor.
 *
 *@params: pointer to struct of parameters *params
 *
 *@returns: the component chi_31, evaluated using the parameters in struct p
 */
double chi_31(struct params * params)
{
  if(params->use_spline == 1)
  {
  //  return chi_13_spline(params);
  }

  return chi_13(params);
}
