
#include "auxiliary/headers.hpp"



int main (int argc, char* argv[]) // command-line input: filename_begin, fileformat, impact_parameter, # of files
{
	typedef double number_type;
	typedef std::size_t size_type;

	// Evaluate command-line input
	std::string filename = argv[1];
	std::string impact_parameter = argv[3];
	filename += impact_parameter;
	filename += "/";
	std::string fileformat = argv[2];
	size_type n_files = to_size_t(argv[4]);


	// Read in and pre-process Trento data
	Collision<number_type> PbPb(10, .2); // Create Collision object

	PbPb.read_in(filename, fileformat, n_files); // read in Trento event files

	PbPb.normalize(1); // normalize events so that integral = 1


	PbPb.centralize(); // shift data so that barycentre at origin

	// set interpolation method: gsl_interp2d_bicubic or gsl_interp2d_bilinear
	const gsl_interp2d_type* interpolation_method = gsl_interp2d_bicubic;
	PbPb.initialize_interpolation(interpolation_method); // initialize objects needed for interpolation


	/*
	goal generate text file with following columns
	m 	l 	e_ml	err_ml
	*/
	// Create outfile name
	std::string outfile = "output/decomposition";
	outfile += impact_parameter;
	outfile += ".txt";
	size_type mMax = 10;
	size_type lMax = 10;
	PbPb.FourierBesselDecompose(outfile, mMax, lMax);

	// // compute phi-averaged energy density profile and save to text file
	// gsl_vector* radii = gsl_vector_alloc(500); // compute energy density for 500 radii between 0 and grid_max-1fm
	// gsl_vector* energy = gsl_vector_alloc(radii->size); // for each radius, the energy density and...
	// gsl_vector* energy_err = gsl_vector_alloc(radii->size); /// ...its statistical uncertainty is saved

	// PbPb.average_azimuthal(radii, energy, energy_err);
	


	// // free allocated memory
	// gsl_vector_free(radii);
	// gsl_vector_free(energy);
	// gsl_vector_free(energy_err);

	PbPb.free();




	return 0;
}