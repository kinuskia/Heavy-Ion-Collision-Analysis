#include <iostream>
#include <cmath>
#include <vector>
#include <gsl/gsl_integration.h>
#include "../auxiliary/to_file.hpp"
#include "glauber_auxiliary.hpp"

int main () 
{
	// Nucleus parameters Pb
	double R = 4.2;
	double a = 0.596;
	double w = 0.;
	double rho0 = get_rho0(R, a, w);
	double N = 208;

	// inelastic cross section
	double sigma_nn_inel = 6.4;
	

	int N_values = 30;
	double b_max = 22;
	std::vector<double> b(N_values);
	std::vector<double> d_sigma_d_b(N_values);

	for (int i = 0; i < N_values; ++i)
	{
		b[i] = b_max * i / N_values;
		d_sigma_d_b[i] = sigma_AB(b[i], N, sigma_nn_inel, rho0, R, a, w);
	}
	std::vector<std::vector<double>> data(2);
	data[0] = b;
	data[1] = d_sigma_d_b;

	to_file("d_sigma_d_b.txt", data);




	return 0;
}