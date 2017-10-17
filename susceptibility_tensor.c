#include <stdio.h>
#include <math.h>
#include "susceptibility_tensor.h"
#include <time.h>
#include <string.h>
#include <stdlib.h>

/*set_params: sets values for the constants (permittivity of free space,
 *            electron charge, electron mass, etc.) as well as the 
 *            particle distribution function and distribution function
 *            parameters. 
 *
 *@params: pointer to struct of parameters *p
 * 
 *@returns: 1 to indicate success //TODO: should we just make this void? 
 */
int set_params(struct params *p)
{
  p->epsilon0  = 1./(4. * M_PI); //permittivity of free space, CGS units
  p->e         = 4.80320680e-10; //electron charge
  p->m         = 9.1093826e-28;  //electron mass
  p->c         = 2.99792458e10;  //speed of light
  p->epsilon   = -1.;            //sign of electron charge
  
  //parameters
  p->B       = 1.;          //background B strength
  p->n_e     = 1.;          //electron number density cm^-3
  p->theta   = M_PI/3.;     //observer angle
  
  //derived quantities
  p->omega_p = sqrt(p->n_e*p->e*p->e / (p->m * p->epsilon0));//plasma frequency    
  p->omega_c = p->e * p->B / (p->m * p->c);               //cyclotron frequency
  
  //integrator parameters
  p->gamma             = 1.5; //will get reset later in integration
  p->real              = 1;   //real part = 1, imag part = 0
  
  //distribution function
  p->dist              = 0; //MJ=0, PL=1, kappa=2
  
  //distribution function parameters
  p->theta_e     = 10.;         //dimensionless electron temp
  p->pl_p        = 3.;          //power-law index, p
  p->gamma_min   = 1.;          //power-law gamma_min
  p->gamma_max   = 1000.;       //power-law gamma_max
  p->kappa       = 3.5;         //kappa index
  p->kappa_width = 10.;         //kappa width, like theta_e
  p->gamma_cutoff = 1e10;       //currently unused

  p->pull_out_Df = 0;
  p->use_spline  = 1;
 
  return 1;
}

/*alpha_I: returns the absorption coefficient alpha_I, for the total intensity
 *         of light along the ray in question, for the given values of 
 *         parameters within the struct p.
 *
 *@params: pointer to struct of parameters *p
 *
 *@returns: absorption coefficient for total intensity (Stokes I) 
 */
double alpha_I(struct params *p)
{
  p->real          = 0;
  double prefactor = 2. * M_PI * p->epsilon0 * p->omega / p->c;
  double term11    = (chi_11(p) * pow(cos(p->theta), 2.)  
  		  + chi_33(p) * pow(sin(p->theta), 2.)
  		  - 2. * chi_13(p) * sin(p->theta) * cos(p->theta));
  double term22    = chi_22(p);
  double ans       = prefactor * (term11 + term22);
  return ans;
}

/*alpha_Q: returns the absorption coefficient alpha_Q, for linearly polarized
 *         light along the ray in question, for the given values of parameters 
 *         within the struct p.
 *
 *@params: pointer to struct of parameters *p
 *
 *@returns: absorption coefficient for linearly polarized light (Stokes Q) 
 */
double alpha_Q(struct params *p)
{
  p->real          = 0;
  double prefactor = 2. * M_PI * p->epsilon0 * p->omega / p->c;
  double term11    = (chi_11(p) * pow(cos(p->theta), 2.)
                    + chi_33(p) * pow(sin(p->theta), 2.)
                    - 2. * chi_13(p) * sin(p->theta) * cos(p->theta));
  double term22    = chi_22(p);
  double ans       = prefactor * (term11 - term22);
  return ans;
}

/*rho_Q: returns the Faraday conversion coefficient rho_Q, which corresponds
 *       to the conversion between linearly polarized and circularly
 *       polarized light by the medium.  The coefficient is calculated for
 *       the given values of parameters within the struct p.
 *
 *@params: pointer to struct of parameters *p
 *
 *@returns: Faraday conversion coefficient rho_Q 
 */
double rho_Q(struct params *p)
{
  p->real          = 1;
  double prefactor = 2. * M_PI * p->epsilon0 * p->omega / p->c;
  double term11    = (chi_11(p) * pow(cos(p->theta), 2.)
                    + chi_33(p) * pow(sin(p->theta), 2.)
                    - 2. * chi_13(p) * sin(p->theta) * cos(p->theta));
  double term22    = chi_22(p);
  double ans       = prefactor * (term22 - term11);
  return ans;
}

/*alpha_V: returns the absorption coefficient alpha_V, for the circularly
 *         polarized light along the ray in question, for the given values of 
 *         parameters within the struct p.  Uses the IEEE/IAU convention for
 *         the sign of Stokes V.
 *
 *@params: pointer to struct of parameters *p
 *
 *@returns: absorption coefficient for circularly polarized light (Stokes V) 
 */
double alpha_V(struct params *p)
{
  p->real            = 1;
  double prefactor   = 4. * M_PI * p->epsilon0 * p->omega / p->c;
  double term1     = (chi_12(p) * cos(p->theta) - chi_32(p) * sin(p->theta));
  double ans       = prefactor * term1;
  return ans;
}

/*rho_V: returns the Faraday rotation coefficient rho_V, which rotates the
 *       plane of polarization (EVPA) for linearly polarized light,
 *       for the given values of parameters within the struct p.
 *
 *@params: pointer to struct of parameters *p
 *
 *@returns: Faraday rotation coefficient rho_V 
 */
double rho_V(struct params *p)
{
  p->real          = 0;
  double prefactor = 4. * M_PI * p->epsilon0 * p->omega / p->c;
  double term1     = (chi_12(p) * cos(p->theta) - chi_32(p) * sin(p->theta));
  double ans       = prefactor * term1;
  return ans;
}

/*plotter: prints the values of the gamma integrand for the component of chi_ij
 *         determined by p.tau_integrand, from gamma=start to gamma=end, in
 *         increments of step.  These values are printed to a file called
 *         output.txt, and can be plotted easily by an external plotting
 *         software to determine if the gamma integrand is being properly
 *         resolved.
 *
 *@params: struct of parameters p
 *
 *@returns: 0 when completed and prints the gamma integrand to a file for
 *          plotting //TODO: make this function return nothing? 
 */
double plotter(struct params p)
{
  FILE *fp;
  fp = fopen("output.txt", "w");
  
  double start = 1.;
  double end   = 1.01;
  double i     = start;
  double step  = 0.00001;
  
  p.tau_integrand = &chi_12_integrand;
  
  while(i < end)
  {
    fprintf(fp, "\n%e    %e", i, gamma_integrand(i, &p));
    printf("\n%e", i);
    i = i + step;
  }
  printf("\n");
  
  return 0.;
}

double spline_plotter(struct params p)
{

  p.pull_out_Df = 1;

  int step_array_size = 271;

  char ch, buffer[1000];
  int i = 0, j = 0;
  float arr[step_array_size];  

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

//  for(i = 0; i < step_array_size; i++)
//  {
//    printf("\n%f", arr[i]);
//  }

  
  FILE *fp;
  fp = fopen("chi_12_real_mod_step.txt", "w");

 
  double start = 1.;
  i        = 0;
  j        = 0;
  double gamma = 1.;

  int max_index = step_array_size;

  p.tau_integrand = &chi_12_integrand;
  
  int array_width = max_index;

  double gamma_omratio_array[array_width][array_width];

  printf("\narray is %d^2 in size\n", array_width);

  for(i = 0; i < max_index; i++)
  {
    p.omega = arr[i] * p.omega_c;

    #pragma omp parallel for
    for(j = 0; j < max_index; j++)
    {
      gamma = arr[j];
      gamma_omratio_array[i][j] = gamma_integrand(gamma, &p);
    }
      printf("\nrow: %d", i);
  }
  printf("\n");

  for(i = 0; i < max_index; i++)
  {
    for(j = 0; j < max_index; j++)
    {
      fprintf(fp, "	%e", gamma_omratio_array[i][j]);
    }
    fprintf(fp, "\n");
  }
  
  return 0.;
}

/*main: sets parameters, runs some calculation, and prints the CPU time elapsed
 *
 *@params: none
 *
 *@returns: nothing
 */
int main(void)
{
  struct params p;
  
  /*set parameters*/
  set_params(&p);
  p.omega = 1. * p.omega_c;
  p.real  = 1;
  
  /*print gamma	gamma_integrand(gamma) with the function plotter(params)*/
//  spline_plotter(p);
  
  /*print omega/omega_c	alpha_S(params)*/
//  printf("\n%e    %e\n", p.omega/p.omega_c, alpha_V(&p));
//  p.tau_integrand = &chi_12_integrand;
//  printf("\n%e    %e\n", p.omega/p.omega_c, gamma_integrand(1.1, &p));
}
