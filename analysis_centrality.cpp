
#include "auxiliary/headers.hpp"



int main (int argc, char* argv[]) // command-line input: filename_begin, fileformat,  # of files
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

	PbPb.read_in(filename, fileformat, n_files, true); // read in Trento event files

	PbPb.normalize(1); // normalize events so that integral = 1


	PbPb.centralize(); // shift data so that barycentre at origin

	//set interpolation method: gsl_interp2d_bicubic or gsl_interp2d_bilinear
	const gsl_interp2d_type* xy_interpolation_method = gsl_interp2d_bicubic;
	PbPb.initialize_xy_interpolation(xy_interpolation_method); // initialize objects needed for interpolation

	std::time_t current_time = std::time(nullptr);
	
	std::cout << current_time-start << "s: " << "Data has been read in and normalized.\n"; 

	
	// Print file with the following columns: impact parameter, number of participants, multiplicity
	PbPb.collision_specs_to_file("output/collision_specs.txt");

	// compute the respective multiplicity limits for specific centrality classes
	std::vector<number_type> classes(12);
	classes[0] = 0;
	classes[1] = 5;
	classes[2] = 10;
	classes[3] = 20; // this generates the classes 0-5%, 5-10, 10-20, 20-30, ..., 80-90, 90-100
	classes[4] = 30;
	classes[5] = 40;
	classes[6] = 50;
	classes[7] = 60;
	classes[8] = 70;
	classes[9] = 80;
	classes[10] = 90;
	classes[11] = 100;
	PbPb.get_percentiles(classes);


	//compute two-point correlation functions for each centrality class

	const gsl_interp_type* r_interpolation_method = gsl_interp_cspline;
	PbPb.initialize_two_point_evaluations(r_interpolation_method, start);

	for (int m = 0; m < 10; ++m)
	{
		for (int c = 1; c < classes.size(); ++c)
		{
			std::cout << "m: " << m << " class: " << classes[c] << "%" << "\n";
			// Create outfile name
			std::string centrality_class = std::to_string(int(classes[c-1])) + "-" + std::to_string(int(classes[c]));
			std::string outfile = "output/two_point_" + centrality_class;
			outfile += "_m";
			outfile += std::to_string(m);
			outfile += ".txt";
			size_type lMax = 20;

			complex_matrix<number_type> TwoPointFunction(lMax, lMax);
			complex_matrix<number_type> TwoPointFunction_err(lMax, lMax);
			PbPb.getTwoPointFunction(m, c, TwoPointFunction, TwoPointFunction_err, r_interpolation_method, start);


			// save moduli of coeffs in text file 
			gsl_matrix* moduli = gsl_matrix_alloc(TwoPointFunction.rowsize(), TwoPointFunction.colsize());
			for (size_type i = 0; i < moduli->size1; ++i)
			{
				for (size_type j = 0; j < moduli->size2; ++j)
				{
					number_type real = TwoPointFunction.get_real(i, j);
					number_type imag = TwoPointFunction.get_imag(i, j);
					gsl_matrix_set(moduli, i, j, sqrt(real*real+imag*imag));
				}
			} 
			to_file(outfile, moduli);
			gsl_matrix_free(moduli);

		}
	}
	

	


	PbPb.free();

	current_time = std::time(nullptr);
	std::cout << current_time-start << "s: " << "Output data saved.\n"; 




	return 0;
}