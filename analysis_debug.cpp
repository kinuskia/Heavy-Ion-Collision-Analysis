
#include "auxiliary/headers.hpp"



int main (int argc, char* argv[]) // command-line input: filename_begin, fileformat, # of files
{
	typedef double number_type;
	typedef std::size_t size_type;

	// number_type dummy;
	// gsl_matrix* profile = gsl_matrix_alloc(100, 100);
	// read_data("Trento/PbPb10000/1000.dat", profile, dummy, dummy, dummy); 

	// for (size_type i = 0; i < profile->size1; ++i)
	// {
	// 	for (size_type j = 0; j < profile->size2; ++j)
	// 	{
	// 		std::cout << gsl_matrix_get(profile, i, j) << " ";
	// 	}
	// 	std::cout << "\n";
	// }

	// Collision<number_type> PbPb(10, 0.2);
	// PbPb.read_in("Trento/PbPb1000/", ".dat", 10000, 4./3); 

	// for (size_type k = 0; k < 10000; ++k)
	// {
	// 	std::cout << PbPb.get_filename("", ".dat", k, 10000) << "\n";
	// }






	return 0;
}