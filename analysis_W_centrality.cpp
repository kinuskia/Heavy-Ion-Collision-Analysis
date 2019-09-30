
#include "auxiliary/headers.hpp"



int main (int argc, char* argv[]) // command-line input: filename_begin, fileformat, # of files
{
	typedef double number_type;
	typedef std::size_t size_type;

	// set timer
	std::time_t start = std::time(nullptr);

	// Evaluate command-line input
	std::string filename = argv[1];
	filename += "/";
	std::string fileformat = argv[2];
	size_type n_files = to_size_t(argv[3]);


	// Read in and pre-process Trento data
	Collision<number_type> PbPb(10, .2); // Create Collision object

	PbPb.read_in(filename, fileformat, n_files, 4./3); // read in Trento event files


	PbPb.centralize(); // shift data so that barycentre at origin


	// define the respective multiplicity limits for specific centrality classes
	std::vector<number_type> classes(101);
	for (size_type i = 0; i < classes.size(); ++i)
	{
		classes[i] = i; // centrality classes from 0 to 100 in 1% intervals
	}
	PbPb.get_percentiles(classes);

	std::time_t current_time = std::time(nullptr);
	std::cout << current_time-start << "s: " << "Data has been read in. \n"; 



	// set interpolation method: gsl_interp2d_bicubic or gsl_interp2d_bilinear
	const gsl_interp2d_type* xy_interpolation_method = gsl_interp2d_bicubic;
	PbPb.initialize_xy_interpolation(xy_interpolation_method); // initialize objects needed for interpolation

	
	
	//PbPb.initialize_n_point_evaluations(r_interpolation_method, start);
	/* Fourier-decompose the profiles */
	PbPb.decompose_azimuthal();
	PbPb.getNormalizations();

	// print weighting functions W
	PbPb.print_W("output/weight_functions", false);

	// Print file with the following columns: impact parameter, number of participants, multiplicity
	PbPb.collision_specs_to_file("output/collision_specs.txt");






	return 0;
}