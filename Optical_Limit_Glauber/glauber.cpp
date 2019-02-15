#include <iostream>
#include <cmath>
#include <vector>
#include <gsl/gsl_integration.h>
#include "../auxiliary/to_file.hpp"
#include "glauber_auxiliary.hpp"

int main (int argc, char* argv[]) // input: x, filename
{
	// Nucleus parameters Pb
	double R = 6.62;
	double a = 0.546;
	double w = 0.;
	double rho0 = get_rho0(R, a, w);
	
	// collision parameters
	double x = std::stod(argv[1])/100.; //  x= 0: wounded nucleon model
	double mult_nn = 1; // average nucleon-nucleon multiplicity
	int N = 207; // collision nucleon number
	double sigma_nn_inel = 6.4; // inelastic nucleon-nucleon cross section in fm^2

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

	//to_file("mult-b_x100.txt", data);
	std::string filename = argv[2]; 
	filename += ".txt";
	to_file(filename, data);




	return 0;
}