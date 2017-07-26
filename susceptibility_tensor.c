#include <stdio.h>
#include <math.h>
#include "susceptibility_tensor.h"
#include <time.h>

int set_params(struct params *p)
{
	p->epsilon0  = 1./(4. * M_PI);
	p->e         = 4.80320680e-10;
	p->m         = 9.1093826e-28;
	p->c         = 2.99792458e10;
	p->epsilon   = -1.;        //sign of electron charge
	
	//parameters
	p->B       = 1.;          //background B strength
	p->n_e     = 1.;          //electron number density cm^-3
	p->theta   = M_PI/3.;     //observer angle

	//derived quantities
	p->omega_p = sqrt(p->n_e * p->e*p->e / (p->m * p->epsilon0));       //plasma frequency    
        p->omega_c = p->e * p->B / (p->m * p->c);                        //cyclotron frequency

	//integrator parameters
	p->gamma             = 1.5; //will get reset later in integration
	p->real              = 1;

	//distribution function
	p->dist              = 0;

	//distribution function parameters
	p->theta_e   = 10;         //dimensionless electron temp
	p->pl_p      = 3;           //power-law index, p
	p->gamma_min = 1.;          //power-law gamma_min
	p->gamma_max = 1000.;       //power-law gamma_max

	return 1;
}

double alpha_I(struct params *p)
{
	p->real          = 0;
        double prefactor = 2. * M_PI * p->epsilon0 * p->omega / p->c;
        double term11    = (chi_11(p) * pow(cos(p->theta), 2.)  
			  + chi_33(p) * pow(sin(p->theta), 2.)
			  + 2. * chi_13(p) * sin(p->theta) * cos(p->theta));
	double term22    = chi_22(p);
        double ans       = prefactor * (term11 + term22);
        return ans;
}

double alpha_Q(struct params *p)
{
        p->real          = 0;
        double prefactor = 2. * M_PI * p->epsilon0 * p->omega / p->c;
        double term11    = (chi_11(p) * pow(cos(p->theta), 2.)
                          + chi_33(p) * pow(sin(p->theta), 2.)
                          + 2. * chi_13(p) * sin(p->theta) * cos(p->theta));
        double term22    = chi_22(p);
        double ans       = prefactor * (term11 - term22);
        return ans;
}

double rho_Q(struct params *p)
{
        p->real          = 1;
        double prefactor = 2. * M_PI * p->epsilon0 * p->omega / p->c;
        double term11    = (chi_11(p) * pow(cos(p->theta), 2.)
                          + chi_33(p) * pow(sin(p->theta), 2.)
                          + 2. * chi_13(p) * sin(p->theta) * cos(p->theta));
        double term22    = chi_22(p);
        double ans       = prefactor * (term22 - term11);
        return ans;
}

double alpha_V(struct params *p)
{
	p->real          = 1;
	double prefactor = 4. * M_PI * p->epsilon0 * p->omega / p->c;
	double term1     = (chi_12(p) * cos(p->theta) + chi_32(p) * sin(p->theta));
	double ans       = prefactor * term1;
	return ans;
}

double rho_V(struct params *p)
{
        p->real          = 0;
        double prefactor = 4. * M_PI * p->epsilon0 * p->omega / p->c;
        double term1     = (chi_12(p) * cos(p->theta) + chi_32(p) * sin(p->theta));
        double ans       = prefactor * term1;
        return ans;
}

double plotter(struct params p)
{
	FILE *fp;
	fp = fopen("output.txt", "w");

	double start = 1.;
	double end   = 150.;
	double i     = start;
	double step  = 2.;

	p.integrand = chi_11_integrand;

        while(i < end)
        {
                fprintf(fp, "\n%e    %e", i, tau_integrator(i, &p));
		printf("\n%e", i);
                i = i + step;
        }
        printf("\n");

	return 0.;
}

int main(void)
{
	/*start timer*/
	clock_t start = clock(), diff;
        struct params p;

	/*set parameters*/
	set_params(&p);
	p.omega = 10000. * p.omega_c;
	p.real  = 0;

	/*print gamma	gamma_integrand(gamma) with the function plotter(params)*/
	plotter(p);

	/*print omega/omega_c	alpha_I(params)*/
//	printf("\n%e    %e\n", p.omega/p.omega_c, alpha_I(&p));
//	printf("\n%e    %e\n", p.omega/p.omega_c, chi_12(&p));	

	/*calculate and print elapsed time*/
	diff = clock() - start;
	int msec = diff * 1000 / CLOCKS_PER_SEC;
	printf("Time taken %d seconds %d milliseconds\n", msec/1000, msec%1000);
}

