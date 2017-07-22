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
        int resolution_factor;
	int real;
	int dist;
};

double I_1_of_2(double alpha, double delta);
double I_1_analytic(double alpha, double delta);
double I_2_analytic(double alpha, double delta);
double MJ(struct params * params);
double Df(struct params * params);

double tau_integrator_11(double gamma, void * parameters);
double tau_integrator_12(double gamma, void * parameters);
double tau_integrator_13(double gamma, void * parameters);
double tau_integrator_32(double gamma, void * parameters);
double tau_integrator_33(double gamma, void * parameters);
double tau_integrator_22(double gamma, void * parameters);

double chi_33_integrand(double tau_prime, void * parameters);

double chi_11(struct params * p);
double chi_22(struct params * p);
double chi_33(struct params * p);
double chi_12(struct params * p);
double chi_21(struct params * p);
double chi_23(struct params * p);
double chi_32(struct params * p);
double chi_13(struct params * p);
double chi_31(struct params * p);
