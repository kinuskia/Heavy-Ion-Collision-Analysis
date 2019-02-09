#include <iostream>
#include <cmath>
#include <vector>
#include <gsl/gsl_integration.h>
#include "../auxiliary/to_file.hpp"
#include "glauber_auxiliary.hpp"

int main ()
{
	// Nucleus parameters
	double R = 6.38;
	double a = 0.535;
	double w = 0.;
	double rho0 = get_rho0(R, a, w);
	
	// collision parameters
	double x = 1.0; //  x= 0: wounded nucleon model
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

	to_file("mult-b_x100.txt", data);





	return 0;
}