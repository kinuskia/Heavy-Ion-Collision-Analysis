#include <iostream>
#include <vector>
//#include <gsl/gsl_vector.h>

#include "model.hpp"
#include "../auxiliary/fbdecomposition_s.hpp" // load after (!) model.hpp
#include "../auxiliary/to_file.hpp"
#include "../auxiliary/to_size_t.hpp"

#include <string>


/*
Idea: In the file "model.hpp" one specifies the position space one-point
and two-point correlation function of an arbitrary initial-state model.
*/

int main (int argc, char* argv[]) // command-line input: centrality_min, centrality_max, IR regulator m, saturation scale Qs0
{
	typedef std::size_t size_type;
	typedef double number_type;

	// Set up initial-state model
	size_type centrality_min = to_size_t(argv[1]);
	size_type centrality_max = to_size_t(argv[2]);
	number_type m_IR = std::stod(argv[3]);
	number_type Qs0 = std::stod(argv[4]);
	Model<number_type> model(m_IR, Qs0);
	std::string centrality = std::to_string(centrality_min) +  "-" + std::to_string(centrality_max);
	model.initialize_W("weight_functions_"+centrality+".txt");


	// Set up Fourier-Bessel decomposition object
	// with rMax = 10 as maximal radial integration length
	FBDecompositionSimplified<number_type> decomposition(model, 9.604);

	decomposition.initialize();

	// Compute <e_l1^(m)e_l2^(-m)> as a function of l
	int mMax = 4;
	int lMax = 10;

	for (int m = mMax; m >= 0; --m)
	{
		// save result in matrix
		gsl_matrix* result = gsl_matrix_alloc(lMax, lMax);
		for (int l1 = 1; l1 <= lMax; ++l1)
		{
			for (int l2 = 1; l2 <= lMax; ++l2)
			{
				std::cout << "m=" << m << ", l1= " << l1 << ", l2=" << l2 << "\n";
			
				number_type current = decomposition.TwoMode(m, l1, -m, l2);	
				gsl_matrix_set(result, l1-1, l2-1, current);
			}
			
		}

		// save result to text file
		std::string filename = "output/"+ centrality +"/two_point_random_connected_m_";
		filename += std::to_string(m);
		//filename += "_test";
		filename += ".txt";

		to_file(filename, result);

		gsl_matrix_free(result);
	}



	return 0;
}