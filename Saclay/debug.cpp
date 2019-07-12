#include <iostream>

#include "model.hpp"

#include <vector>
#include<cmath>
#include "../auxiliary/to_file.hpp"

#include <gsl/gsl_sf_bessel.h>

int main ()
{
	typedef std::size_t size_type;
	typedef double number_type;
	Model<number_type> model(1e-1);
	model.initialize_W("weight_functions.txt");

	size_type N = 100;
	std::vector<number_type> x(N);
	std::vector<number_type> f(N);
	std::vector<number_type> g(N);

	for (size_type i = 0; i < N; ++i)
	{
		number_type r_min = 0;
		number_type r_max = 9;
		number_type phi1 = 3.927;
		number_type phi2 = 6.234;
		number_type r2 = 8.022;
		x[i] = r_min + (r_max- r_min)*i/N;
		f[i] = model.TwoPoint(x[i]*cos(phi1), x[i]*sin(phi1), r2*cos(phi2), r2*sin(phi2));
		g[i] = model.OnePoint(x[i]*cos(phi1), x[i]*sin(phi1));
	}

	std::vector<std::vector<number_type>> data(3);
	data[0] = x;
	data[1] = f;
	data[2] = g;

	to_file("debug.txt", data);

	// size_type N = 100;
	// std::vector<number_type> x(N);
	// std::vector<number_type> f(N);
	// std::vector<number_type> g(N);

	// for (size_type i = 0; i < N; ++i)
	// {
	// 	number_type r_min = 0.1;
	// 	number_type r_max = 20;
	// 	number_type m = 1e-1;
	// 	x[i] = r_min + (r_max- r_min)*i/N;
	// 	f[i] = gsl_sf_bessel_K1(m*x[i]);
	// 	g[i] = log(4./m/m/x[i]/x[i]);
	// }

	// std::vector<std::vector<number_type>> data(3);
	// data[0] = x;
	// data[1] = f;
	// data[2] = g;

	// to_file("debug.txt", data);




	return 0;
}