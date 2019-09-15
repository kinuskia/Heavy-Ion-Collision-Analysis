#include <iostream>

#include "model.hpp"


#include <vector>
#include <string>
#include "../auxiliary/to_file.hpp"
#include "../auxiliary/fbdecomposition.hpp"

int main ()
{
	typedef double number_type;
	typedef std::size_t size_type;
	
	// Model<number_type> LargeNc(0.14);
	


	// size_type centrality_min = size_type(20);
	// size_type centrality_max = size_type(21);
	// std::string centrality = std::to_string(centrality_min) + "-" + std::to_string(centrality_max);
	// LargeNc.initialize_W("weight_functions_"+centrality+".txt");

	// FBDecomposition<number_type> decomposition(LargeNc, 9.604+0*1.*8.604);
	// decomposition.initialize();

	// // set number of radial grid points per dimension
	// decomposition.set_N_discret(29);

	// // set # gridpoints for FFT
	// decomposition.set_N_discret(64);

	// decomposition.get_Bessel_deriv_zeros(10, 10);

	// number_type rmax = 9.604-0.2;
	// size_type N = 100;
	// int m = 1;
	// int l = 9;

	// gsl_matrix* result = gsl_matrix_alloc(N, N);

	// for (size_type i = 0; i < N; ++i)
	// {
	// 	for (size_type j = 0; j < N; ++j)
	// 	{
	// 		number_type r1 = rmax*i/(N-1);
	// 		number_type r2 = rmax*j/(N-1);
	// 		gsl_matrix_set(result, i, j, decomposition.TwoPoint_mode(m, r1, r2));
	// 	}
	// }

	// to_file("TwoPoint_mode-"+ std::to_string(m) +".txt", result);

	// gsl_matrix* FB = gsl_matrix_alloc(N, N);

	// for (size_type i = 0; i < N; ++i)
	// {
	// 	for (size_type j = 0; j < N; ++j)
	// 	{
	// 		number_type r1 = rmax*i/(N-1);
	// 		number_type r2 = rmax*j/(N-1);
	// 		gsl_matrix_set(FB, i, j, decomposition.weight(m, l, r1)*decomposition.weight(m, l, r2));
	// 	}
	// }

	// to_file("FB_weight_ml_"+std::to_string(m)+"-"+std::to_string(l)+".txt", FB);
 
	// std::cout << decomposition.integ_r_test(300) << "\n";


	size_type N = 5000;
	number_type r_min = 0;
	number_type r_max = 5.e-2;
	size_type centrality_min = size_type(20);
	size_type centrality_max = size_type(21);
	number_type R = 0;

	

	std::vector<number_type> r(N);
	std::vector<number_type> y(N);

	Model<number_type> LargeNc(1.e-3);
	std::string centrality = std::to_string(centrality_min) + "-" + std::to_string(centrality_max);

	LargeNc.initialize_W("weight_functions_"+centrality+".txt");

	
	for (size_type i = 0; i < N; ++i)
	{
		r[i] = r_min + (r_max-r_min)*i/N;
		y[i] = LargeNc.TwoPoint(R, r[i]);
	}

	std::vector<std::vector<number_type>> output(2);
	output[0] = r;
	output[1] = y;

	to_file("TwoPoint_"+centrality+".txt", output);



	return 0;
}