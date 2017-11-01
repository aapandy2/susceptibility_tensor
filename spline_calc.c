#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>

#include "susceptibility_tensor.h"

double chi_11_spline(struct params * params)
{
  params->component = 11;
  return chi_ij_spline(params);
}

double chi_12_spline(struct params * params)
{
  params->component = 12;
  return chi_ij_spline(params);
}

double chi_13_spline(struct params * params)
{
  params->component = 13;
  return chi_ij_spline(params);
}

double chi_22_spline(struct params * params)
{
  params->component = 22;
  return chi_ij_spline(params);
}

double chi_32_spline(struct params * params)
{
  params->component = 32;
  return chi_ij_spline(params);
}

double chi_33_spline(struct params * params)
{
  params->component = 33;
  return chi_ij_spline(params);
}

double chi_ij_spline(struct params * params)
{
  double prefactor = 2. * M_PI * params->omega_p*params->omega_p 
                     / (params->omega * params->omega);
  double start = 1.;
  double end   = end_approx(params);

  double ans = prefactor * gauss_legendre_spline(start, end, params);

  return ans;
}

double chi_ij_spline_integrand(double gamma, struct params * params)
{
  params->gamma = gamma;
  double beta = sqrt(1. - pow(gamma, -2.));

  double dist = Df(params);
  double gam_term = dist * pow(gamma, 3.) * pow(beta, 3.);
  double ans = gam_term * spline_integrand(gamma, 
                                           params->omega/params->omega_c, 
                                           params);

  return ans;
}

#define NUM_QUAD 21
double gauss_legendre_spline(double start, double end, struct params * params)
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

//  # pragma omp parallel for private(i , x, weight) reduction ( + : sum )
  for(i = 0; i < n; i++)
  {
    x        = quadPts[i];
    weight   = weights[i];

    sum = sum + (end - start)/2. 
                * chi_ij_spline_integrand((end - start)/2. * x + (end + start)/2., params) * weight;
  }

  return sum;
}

#define STEP_ARRAY_SIZE 271
//#define STEP_ARRAY_SIZE 200
double spline_integrand(double gamma, double omratio, struct params * params)
{
  static int loaded_file = 0;

  //TODO: do this for all possible changeable parameters
  static int past_component = 0;

  if(past_component == 0 ||
     past_component != params->component)
  {
    loaded_file = 0;
    past_component = params->component;
  }

  char ch, buffer[1000];
  int i = 0, j = 0;
  static double arr[STEP_ARRAY_SIZE];

  double myvariable;

  static double z_vals[STEP_ARRAY_SIZE * STEP_ARRAY_SIZE];

//  char fname[34]; //= "chi_12_real_step_n.txt";
  char fname[36];

  /*choose file*/
  if(params->component == 11)
  {
    if(params->real == 1)
    {
      strcpy(fname, "./datafiles/chi_11_real_mod_step.txt");
    }
    else
    {
      strcpy(fname, "./datafiles/chi_11_imag_mod_step.txt");
    }
  }
  else if(params->component == 12)
  {
    if(params->real == 1)
    {
      strcpy(fname, "./datafiles/chi_12_real_mod_step.txt");
    }
    else
    {
      strcpy(fname, "./datafiles/chi_12_imag_mod_step.txt");
    }
  }
  else if(params->component == 13)
  {
    if(params->real == 1)
    {
      strcpy(fname, "./datafiles/chi_13_real_mod_step.txt");
    }
    else
    {
      strcpy(fname, "./datafiles/chi_13_imag_mod_step.txt");
    }
  }
  else if(params->component == 22)
  {
    if(params->real == 1)
    {
      strcpy(fname, "./datafiles/chi_22_real_mod_step.txt");
    }
    else
    {
      strcpy(fname, "./datafiles/chi_22_imag_mod_step.txt");
    }
  }
  else if(params->component == 32)
  {
    if(params->real == 1)
    {
      strcpy(fname, "./datafiles/chi_32_real_mod_step.txt");
    }
    else
    {
      strcpy(fname, "./datafiles/chi_32_imag_mod_step.txt");
    }
  }
  else if(params->component == 33)
  {
    if(params->real == 1)
    {
      strcpy(fname, "./datafiles/chi_33_real_mod_step.txt");
    }
    else
    {
      strcpy(fname, "./datafiles/chi_33_imag_mod_step.txt");
    }
  }
  else
  {
    printf("\ncorresponding file not found\n");
  }

    if(loaded_file == 0)
  {

  FILE *fs; 
 
  // Openning the file with file handler as fs
  fs = fopen("step_array.txt", "r");
  
  // Read the file unless the file encounters an EOF
  while(1)
  {
   // Reads the character where the seeker is currently
   ch = fgetc(fs);
   
   // If EOF is encountered then break out of the while loop
   if(ch == EOF){
   	break;
   }
   
   // If the delimiter is encountered(which can be
   // anything according to your wish) then skip the character
   // and store the last read array of characters in
   // an integer array
   else if(ch == ' '){
   
   	// Converting the content of the buffer into
   	// an array position
   	arr[j] = atof(buffer);
   
   	// Incrementing the array position
   	j++;
   
   	// Clearing the buffer, this function takes two
   	// arguments, one is a character pointer and 
   	// the other one is the size of the character array
   	bzero(buffer, 32);
   
   	// clearing the counter which counts the number
   	// of character in each number used for reading
   	// into the buffer.
   	i = 0;
   
   	// then continue
   	continue;
   }
   else{
   
   	// reads the current character in the buffer
   	buffer[i] = ch;
   
   	// increamenting the counter so that the next
   	// character is read in the next position in 
   	// the array of buffer
   	i++;
   }
    }

//  for(i = 0; i < STEP_ARRAY_SIZE; i++)
//  {
//    printf("\n%f", arr[i]);
//  }

    FILE *myfile;
    myfile=fopen(fname, "r");

    for(i = 0; i < STEP_ARRAY_SIZE; i++)
    {
      for (j = 0 ; j < STEP_ARRAY_SIZE; j++)
      {
        fscanf(myfile,"%lf", &myvariable);
        z_vals[j*STEP_ARRAY_SIZE + i] = myvariable;
      }
    }

    fclose(myfile);

    loaded_file = 1;
  }


//  const gsl_interp2d_type *T = gsl_interp2d_bilinear;
  const gsl_interp2d_type *T = gsl_interp2d_bicubic;
  gsl_spline2d *spline = gsl_spline2d_alloc(T, STEP_ARRAY_SIZE, STEP_ARRAY_SIZE);
  gsl_interp_accel *xacc = gsl_interp_accel_alloc();
  gsl_interp_accel *yacc = gsl_interp_accel_alloc();

  /* initialize interpolation */
  gsl_spline2d_init(spline, arr, arr, z_vals, STEP_ARRAY_SIZE, STEP_ARRAY_SIZE);

  double ans = gsl_spline2d_eval(spline, omratio, gamma, xacc, yacc);

  gsl_spline2d_free(spline);
  gsl_interp_accel_free(xacc);
  gsl_interp_accel_free(yacc);

  return ans;
}
