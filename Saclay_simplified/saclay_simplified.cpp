#include <iostream>
#include <vector>
//#include <gsl/gsl_vector.h>

#include "model.hpp"
#include "../auxiliary/fbdecomposition_s.hpp" // load after (!) model.hpp
#include "../auxiliary/to_file.hpp"



/*
Idea: In the file "model.hpp" one specifies the position space one-point
and two-point correlation function of an arbitrary initial-state model.
*/

int main ()
{
	typedef std::size_t size_type;
	typedef double number_type;

	// Set up initial-state model
	Model<number_type> model(6e-2);

	model.initialize_W("weight_functions.txt");


	// Set up Fourier-Bessel decomposition object
	// with rMax = 10 as maximal radial integration length
	FBDecompositionSimplified<number_type> decomposition(model, 9.604);

	decomposition.initialize();

	// Compute <e_l1^(m)e_l2^(-m)> as a function of l
	int mMax = 4;
	int lMax = 6;

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
		std::string filename = "output/two_point_random_connected_m_";
		filename += std::to_string(m);
		//filename += "_test";
		filename += ".txt";

		to_file(filename, result);

		gsl_matrix_free(result);
	}



	return 0;
}