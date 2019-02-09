
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

	PbPb.read_in(filename, fileformat, n_files); // read in Trento event files

	PbPb.normalize(1); // normalize events so that integral = 1


	PbPb.centralize(); // shift data so that barycentre at origin

	// set interpolation method: gsl_interp2d_bicubic or gsl_interp2d_bilinear
	const gsl_interp2d_type* xy_interpolation_method = gsl_interp2d_bicubic;
	PbPb.initialize_xy_interpolation(xy_interpolation_method); // initialize objects needed for interpolation

	std::time_t current_time = std::time(nullptr);
	std::cout << current_time-start << "s: " << "Data has been read in and normalized.\n"; 




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

	/*
	compute expectation values <e_0l> for each centrality class (m=0 for non-vanishing values)
	*/

	const gsl_interp_type* r_interpolation_method = gsl_interp_cspline;
	PbPb.initialize_n_point_evaluations(r_interpolation_method, start);

	// print weighting functions W
	PbPb.print_W("output/weight_functions.txt", 200);

	for (int c = 1; c < classes.size(); ++c)
	{
		std::cout << " class: " << classes[c] << "%" << "\n";
		// Create outfile names
		std::string centrality_class = std::to_string(int(classes[c-1])) + "-" + std::to_string(int(classes[c]));
		std::string outfile_modulus = "output/one_point_" + centrality_class;
		outfile_modulus += "_modulus";
		outfile_modulus += ".txt";
		size_type lMax = 10;

		complex_matrix<number_type> OnePointFunction(1, lMax);
		complex_matrix<number_type> OnePointFunction_err(1, lMax);

		PbPb.getOnePointFunction(c, OnePointFunction, OnePointFunction_err, start);
		
		// save moduli of coeffs in text file 
		std::vector<number_type> moduli(lMax);
		for (size_type i = 0; i < moduli.size(); ++i)
		{
			number_type real = OnePointFunction.get_real(0, i);
			number_type imag = OnePointFunction.get_imag(0, i);
			moduli[i] = sqrt(real*real+imag*imag);	
		}

		std::vector<std::vector<number_type>> moduli_column(1);
		moduli_column[0] = moduli; 
		to_file(outfile_modulus, moduli_column);

		// save phase of coeffs in text file
		std::string outfile_phase = "output/one_point_" + centrality_class;
		outfile_phase += "_phase";
		outfile_phase += ".txt";
		std::vector<number_type> phases(lMax);
		for (size_type i = 0; i < phases.size(); ++i)
		{
			number_type real = OnePointFunction.get_real(0, i);
			number_type imag = OnePointFunction.get_imag(0, i);
			phases[i] = atan2(imag, real);	
		}
		std::vector<std::vector<number_type>> phases_column(1);
		phases_column[0] = phases;

		to_file(outfile_phase, phases_column); 

	}




	current_time = std::time(nullptr);
	std::cout << current_time-start << "s: " << "Output data saved.\n"; 




	return 0;
}