struct params
{
        double epsilon0;
	double epsilon;
	double e;
	double m;
	double c;
	double B;
	double n_e;
	double theta;
	double theta_e;
	double pl_p;
	double gamma_min;
	double gamma_max;
	double omega_c;
	double omega_p;
        double omega;
        double gamma;
	int real;
	int dist;
	double (*integrand)(double, void * parameters);
};

double I_1_of_2(double alpha, double delta);
double I_1_analytic(double alpha, double delta);
double I_2_analytic(double alpha, double delta);
double MJ(struct params * params);
double Df(struct params * params);

double chi_11_integrand(double tau_prime, void * parameters);
double chi_12_integrand(double tau_prime, void * parameters);
double chi_32_integrand(double tau_prime, void * parameters);
double chi_13_integrand(double tau_prime, void * parameters);
double chi_33_integrand(double tau_prime, void * parameters);
double chi_22_integrand_p1(double tau_prime, void * parameters);
double chi_22_integrand_p2(double tau_prime, void * parameters);
double chi_22_integrand_real(double tau_prime, void * parameters);

double tau_integrator(double gamma, void * parameters);
double gamma_integrator(struct params * p);

double chi_11(struct params * p);
double chi_22(struct params * p);
double chi_33(struct params * p);
double chi_12(struct params * p);
double chi_21(struct params * p);
double chi_23(struct params * p);
double chi_32(struct params * p);
double chi_13(struct params * p);
double chi_31(struct params * p);
