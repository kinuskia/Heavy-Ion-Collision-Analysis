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
	

	int N_values = 30;
	double r_max = 12;
	std::vector<double> r(N_values);
	std::vector<double> density(N_values);

	for (int i = 0; i < N_values; ++i)
	{
		r[i] = r_max * i / N_values;
		density[i] = rho(r[i], rho0, R, a, w);
	}
	std::vector<std::vector<double>> data(2);
	data[0] = r;
	data[1] = density;

	to_file("Cu_density.txt", data);




	return 0;
}