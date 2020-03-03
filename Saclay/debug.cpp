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
	
	// Model<number_type> LargeNc(1e-3);
	


	// size_type centrality_min = size_type(20);
	// size_type centrality_max = size_type(21);
	// std::string centrality = std::to_string(centrality_min) + "-" + std::to_string(centrality_max);
	// LargeNc.initialize_W("weight_functions_"+centrality+".txt");

	// FBDecomposition<number_type> decomposition(LargeNc, 9.604+0*1.*8.604);
	// decomposition.initialize();

	// // set number of radial grid points per dimension
	// decomposition.set_N_discret(41);

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
 
 // 	size_type Nmax = 100;
 // 	for (size_type i = 0; i < Nmax; ++i)
 // 	{
 // 		number_type r = (9.604-0.2)*(i+1)/Nmax;
 // 		number_type y = decomposition.weight(1, 9, r)*decomposition.weight(1, 9, r)*decomposition.TwoPoint_mode(1, r, r);
 // 		std::cout << r << " " << y << "\n";
 // 	}




size_type N = 5000;
	number_type r_min = 0.1;
	number_type r_max = 9;
	size_type centrality_min = size_type(0);
	size_type centrality_max = size_type(1);
	number_type R = 0;

	

	std::vector<number_type> r(N);
	std::vector<number_type> y(N);

	Model<number_type> LargeNc(0.04);
	std::string centrality = std::to_string(centrality_min) + "-" + std::to_string(centrality_max);

	LargeNc.initialize_W("weight_functions_"+centrality+".txt");

	for (size_type i = 0; i < 100; ++i)
	{
		number_type r = r_min + (r_max-r_min)*i/100;
		number_type alpha_s = 0.4095*0+0.25;
		number_type g = sqrt(4.*3.1415925*alpha_s);
		number_type coupling = g;
		number_type Q0 = 1.24;
		number_type W0 = 0.040649297650144;
		number_type IR = 0.14;
		number_type hc = 0.1973;
		//number_type Qs2 = LargeNc.Q2(r, 1.24, 0.040649297650144)*LargeNc.modifiedGamma(0.14, r/hc) *8.*3.1415926;
		//std::cout << r << " " << LargeNc.modifiedGamma(0.14, r/hc) << "\n";
		//std::cout << r << " " << LargeNc.TwoPoint(r, r) << "\n";
	}
	//std::cout << LargeNc.TwoPoint << "\n";
	// for (size_type i = 0; i < N; ++i)
	// {
	// 	r[i] = r_min + (r_max-r_min)*i/N;
	// 	y[i] = LargeNc.TwoPoint(R, r[i]);
	// }

	// std::vector<std::vector<number_type>> output(2);
	// output[0] = r;
	// output[1] = y;

	// to_file("TwoPoint_"+centrality+".txt", output);



	return 0;
}