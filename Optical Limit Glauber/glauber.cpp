#include <iostream>
#include <cmath>
#include <vector>
#include <gsl/gsl_integration.h>
#include "../auxiliary/to_file.hpp"


struct rho_params
{
	double rho0;
	double R;
	double a;
	double w;
};


// Woods-Saxon nuclear density function
double rho(double r, double rho0, double R, double a, double w)
{
	return rho0*(1.+w*r*r/R/R)/(1.+exp((r-R)/a));
}

// gsl wrapper
double rho_gsl(double r, void* params)
{
	rho_params* params_ = (rho_params*) params;
	double rho0 = params_->rho0;
	double R = params_->R;
	double a = params_->a;
	double w = params_->w;
	return rho(r, rho0, R, a, w);
}

// density function times 4 pi r^2
double rho_V(double r, double rho0, double R, double a, double w)
{
	return 4.*3.1415926*r*r*rho(r, rho0, R, a, w);
}

double rho_V_gsl(double r, void* params)
{
	rho_params* params_ = (rho_params*) params;
	double rho0 = params_->rho0;
	double R = params_->R;
	double a = params_->a;
	double w = params_->w;
	return rho_V(r, rho0, R, a, w);
}

struct rho_z_params
{
	double r;
	double rho0;
	double R;
	double a;
	double w;
};

// nuclear density function in cylindrical coordinates
double rho_z(double z, double r, double rho0, double R, double a, double w)
{	
	return rho(sqrt(r*r+z*z), rho0, R, a, w);
}

double rho_z_gsl(double z, void* params)
{
	rho_z_params* params_ = (rho_z_params*) params;
	double r = params_->r;
	double rho0 = params_->rho0;
	double R = params_->R;
	double a = params_->a;
	double w = params_->w;
	rho_params rho_params_ = {rho0, R, a, w};
	return rho_z(z, r, rho0, R, a, w);
}

// get normalization constant so that int rho(r) dV = 1
double get_rho0(double R, double a, double w)
{
	gsl_integration_workspace* workspace = gsl_integration_workspace_alloc(1000);
	
	gsl_function F;
	F.function = &rho_V_gsl;

	double rho0 = 1;

	rho_params nucleus = {rho0, R, a, w};
	F.params = &nucleus;
	double total;
	double total_err;
	gsl_integration_qagiu(&F, 0, 0, 1.e-3, 1000, workspace, &total, &total_err);
	gsl_integration_workspace_free(workspace);

	return 1./total;
}

// one-nucleus thickness function
double T(double r, double rho0, double R, double a, double w)
{
	gsl_integration_workspace* workspace = gsl_integration_workspace_alloc(1000);
	
	gsl_function F;
	F.function = &rho_z_gsl;

	rho_z_params nucleus = {r, rho0, R, a, w};
	F.params = &nucleus;
	double result;
	double result_err;
	gsl_integration_qagi(&F, 0, 1.e-3, 1000, workspace, &result, &result_err);
	gsl_integration_workspace_free(workspace);

	return result;
}

// two-nuclei thickness function integrand
double T2(double x, double y, double b, double rho0, double R, double a, double w)
{
	double rB = sqrt(x*x+y*y);
	double rA = sqrt((b-x)*(b-x)+y*y);
	return T(rB, rho0, R, a, w) * T(rA, rho0, R, a, w);
}

struct T2_params
{
	double y;
	double b;
	double rho0;
	double R;
	double a;
	double w;	
};

double T2_gsl(double x, void* params)
{
	T2_params* params_ = (T2_params*) params;
	double y = params_->y;
	double b = params_->b;
	double rho0 = params_->rho0;
	double R = params_->R;
	double a = params_->a;
	double w = params_->w;
	return T2(x, y, b, rho0, R, a, w);
}

// integrate T2 function with respect to x
double T2_x(double y, double b, double rho0, double R, double a, double w)
{
	gsl_integration_workspace* workspace = gsl_integration_workspace_alloc(1000);
	
	gsl_function F;
	F.function = &T2_gsl;

	T2_params nucleus = {y, b, rho0, R, a, w};
	F.params = &nucleus;
	double result;
	double result_err;
	gsl_integration_qagi(&F, 0, 1.e-3, 1000, workspace, &result, &result_err);
	gsl_integration_workspace_free(workspace);
	return result;
}

struct T2x_params
{
	double b;
	double rho0;
	double R;
	double a;
	double w;	
};

double T2x_gsl(double y, void* params)
{
	T2x_params* params_ = (T2x_params*) params;
	double b = params_->b;
	double rho0 = params_->rho0;
	double R = params_->R;
	double a = params_->a;
	double w = params_->w;
	return T2_x(y, b, rho0, R, a, w);
}

// integrate T2_x with respect to y
double T_AB(double b, double rho0, double R, double a, double w)
{
	gsl_integration_workspace* workspace = gsl_integration_workspace_alloc(1000);
	
	gsl_function F;
	F.function = &T2x_gsl;

	T2x_params nucleus = {b, rho0, R, a, w};
	F.params = &nucleus;
	double result;
	double result_err;
	gsl_integration_qagi(&F, 0, 1.e-3, 1000, workspace, &result, &result_err);
	gsl_integration_workspace_free(workspace);
	return result;
}

// nucleon-nucleus (nucleon number N) cross section at impact parameter b and given inelastic nn cross section
double sigma(double b, int N, double sigma_nn_inel, double rho0, double R, double a, double w) 
{
	return 1.-exp(-1.*N*T(b, rho0, R, a, w)*sigma_nn_inel);
}

// nucleus-nucleus (nucleon number N) cross section at impact parameter b and given inelastic nn cross section
double sigma_AB(double b, int N, double sigma_nn_inel, double rho0, double R, double a, double w) 
{
	return 1.-exp(-1.*N*N*T_AB(b, rho0, R, a, w)*sigma_nn_inel);
}

// average number of nucleon collisions at impact parameter
double N_coll(double b, int N, double sigma_nn_inel, double rho0, double R, double a, double w) 
{
	return 1.*N*N*T_AB(b, rho0, R, a, w)*sigma_nn_inel;
}

// average number of participating nucleons at impact parameter
double N_part(double b, int N, double sigma_nn_inel, double rho0, double R, double a, double w) 
{
	return 2.*N*sigma(b, N, sigma_nn_inel, rho0, R, a, w)/sigma_AB(b, N, sigma_nn_inel, rho0, R, a, w);
}

// average multiplicity at impact parameter b, mixture constant x and average nucleon-nucleon multiplicity
double N_AB(double x, double mult_nn, double b, int N, double sigma_nn_inel, double rho0, double R, double a, double w) 
{
	return ((1.-x)/2 * N_part(b, N, sigma_nn_inel, rho0, R, a, w) + x * N_coll(b, N, sigma_nn_inel, rho0, R, a, w))*mult_nn;
}

int main ()
{
	// Nucleus parameters
	double R = 6.38;
	double a = 0.535;
	double w = 0.;
	double rho0 = get_rho0(R, a, w);
	
	// collision parameters
	double x = 0.9; //  x= 0: wounded nucleon model
	double mult_nn = 1; // average nucleon-nucleon multiplicity
	int N = 197; // collision nucleon number
	double sigma_nn_inel = 4.; // inelastic nucleon-nucleon cross section in fm^2

	int N_values = 30;
	double b_max = 20;
	std::vector<double> b(N_values);
	std::vector<double> mult(b);

	for (int i = 0; i < N_values; ++i)
	{
		b[i] = b_max * i / N_values;
		mult[i] = N_AB(x, mult_nn, b[i], N, sigma_nn_inel, rho0, R, a, w);
	}
	std::vector<std::vector<double>> data(2);
	data[0] = b;
	data[1] = mult;

	to_file("mult-b_x90_float.txt", data);





	return 0;
}