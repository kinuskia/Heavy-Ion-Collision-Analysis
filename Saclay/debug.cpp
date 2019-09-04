#include <iostream>

#include "../Saclay_simplified/model.hpp"

#include <vector>
#include <string>
#include "../auxiliary/to_file.hpp"

int main ()
{
	typedef double number_type;
	typedef std::size_t size_type;
	
	Model<number_type> LargeNc(0.14, 1.24);
	

	size_type N = 10;
	number_type R = 9.604/sqrt(2);
	size_type centrality_min = size_type(20);
	size_type centrality_max = size_type(21);
	std::string centrality = std::to_string(centrality_min) + "-" + std::to_string(centrality_max);
	LargeNc.initialize_W("../Saclay_simplified/weight_functions_"+centrality+".txt");

	const gsl_interp2d_type* xy_interpolation_method = gsl_interp2d_bicubic;
	LargeNc.initialize_OnePoint("../output/profiles_averaged_"+centrality+".txt", 100, xy_interpolation_method, 10, 0.2);
	// LargeNc.initialize_rho();
	// LargeNc.get_Bessel_deriv_zeros(20,20);
	// LargeNc.initialize_coeffs("../output/one_point_"+centrality+".txt", 10, 10);

	// Print profile
	gsl_matrix* result = gsl_matrix_alloc(N-2, N-2);

	for (size_type i = 1; i < N-1; ++i)
	{
		for (size_type j = 1; j < N-1; ++j)
		{
			number_type gridmax = 10;
			number_type gridstep = 0.2;
			number_type min = -gridmax +gridstep*0.51 ;
			number_type max = gridmax -gridstep*0.51 ;
			number_type x = min + (max-min)*j/(N-1);
			number_type y = max - (max-min)*i/(N-1); 
			number_type current = LargeNc.OnePoint(x, y);
			gsl_matrix_set(result, i-1, j-1, current);
			//std::cout << "(" << x << "," << y << ") ";
			std::cout << current << " ";
		}
		std::cout << "\n";
	}

	to_file("OnePoint_Profile_ml_10-10.txt", result);
	LargeNc.print_OnePoint(10);



	// cut through profile
	std::vector<number_type> r(N);
	std::vector<number_type> y(N);

	for (size_type i = 0; i < N; ++i)
	{
		r[i] = -1.0*R +2.*R*i/(N-1);
		//y[i] = LargeNc.TwoPoint_test(r[i]*cos(phi), r[i]*sin(phi), 0, 0);
		y[i] = LargeNc.TwoPoint(r[i], 0);
	}

	std::vector<std::vector<number_type>> output(2);
	output[0] = r;
	output[1] = y;

	to_file("OnePoint_cut_ml_10-10.txt", output);

	// size_type N = 500;
	// number_type r_min = 0.0001;
	// number_type r_max = 9;
	// size_type centrality_min = size_type(0);
	// size_type centrality_max = size_type(1);
	// number_type phi = 10*3.1415926/180;

	

	// std::vector<number_type> r(N);
	// std::vector<number_type> y(N);

	// Model<number_type> LargeNc(0.14);
	// std::string centrality = std::to_string(centrality_min) + "-" + std::to_string(centrality_max);

	// LargeNc.initialize_W("weight_functions_"+centrality+".txt");
	// LargeNc.initialize_rho();
	// LargeNc.get_Bessel_deriv_zeros(20,20);
	// LargeNc.initialize_coeffs("one_point_"+centrality+".txt");
	
	// for (size_type i = 0; i < N; ++i)
	// {
	// 	r[i] = r_min + (r_max-r_min)*i/N;
	// 	//y[i] = LargeNc.TwoPoint_test(r[i]*cos(phi), r[i]*sin(phi), 0, 0);
	// 	y[i] = LargeNc.Bessel(1,2,r[i]);
	// }

	// std::vector<std::vector<number_type>> output(2);
	// output[0] = r;
	// output[1] = y;

	// for (size_type m = 0; m < 6; ++m)
	// {
	// 	for (size_type l = 1; l < 11; ++l)
	// 	{
	// 		std::cout << LargeNc.Bessel_zero(m, l) << " ";
	// 	}
	// 	std::cout << "\n";
	// }

	// to_file("rho.txt", output);




	return 0;
}