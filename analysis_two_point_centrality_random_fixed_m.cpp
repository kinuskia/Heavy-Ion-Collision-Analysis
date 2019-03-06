
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

	// Randomize reaction plane angle
	PbPb.random_reaction_plane(true);

	PbPb.read_in(filename, fileformat, n_files, true); // read in Trento event files

	//PbPb.normalize(1); // normalize events so that integral = 1


	PbPb.centralize(); // shift data so that barycentre at origin

	//set interpolation method: gsl_interp2d_bicubic or gsl_interp2d_bilinear
	const gsl_interp2d_type* xy_interpolation_method = gsl_interp2d_bicubic;
	PbPb.initialize_xy_interpolation(xy_interpolation_method); // initialize objects needed for interpolation

	std::time_t current_time = std::time(nullptr);
	
	std::cout << current_time-start << "s: " << "Data has been read in and normalized.\n"; 

	
	// Print file with the following columns: impact parameter, number of participants, multiplicity
	PbPb.collision_specs_to_file("output/collision_specs.txt");

	// compute the respective multiplicity limits for specific centrality classes
	std::vector<number_type> classes(6);
	classes[0] = 0;
	classes[1] = 5;
	classes[2] = 10;
	classes[3] = 20; // this generates the classes 0-5%, 5-10, 10-20, 20-30, ..., 80-90, 90-100
	classes[4] = 30;
	classes[5] = 40;
	// classes[6] = 50;
	// classes[7] = 60;
	// classes[8] = 70;
	// classes[9] = 80;
	// classes[10] = 90;
	// classes[11] = 100;
	PbPb.get_percentiles(classes);


	//compute two-point correlation functions for each centrality class

	const gsl_interp_type* r_interpolation_method = gsl_interp_cspline;
	PbPb.initialize_n_point_evaluations(r_interpolation_method, start);

	// print weighting functions W
	PbPb.print_W("output/weight_functions.txt", 200);


	// compute two-point correlation functions that include background term <e00 eml>
	number_type max_rel_error = -1.0; // variable to keep track of current relative error
	std::cout << "Computing background correlations.\n";
	for (int c = 1; c < classes.size(); ++c)
	{
		// Create filename
		std::string outfile_background = "output/two_point_random_background_";
		outfile_background += std::to_string(int(classes[c-1])) + "-" + std::to_string(int(classes[c]));
		outfile_background += "_real";
		std::string outfile_background_err = outfile_background;
		outfile_background_err += "_error.txt";
		outfile_background += ".txt";

		number_type mMax = 4;
		number_type lMax = 9;

		complex_matrix<number_type> TwoPointFunction(2*mMax+1, lMax);
		complex_matrix<number_type> TwoPointFunction_err(2*mMax+1, lMax);
		PbPb.getTwoPointFunction_background(c, TwoPointFunction, TwoPointFunction_err, start);

		// save to output
		gsl_matrix* result = gsl_matrix_alloc(TwoPointFunction.rowsize(), TwoPointFunction.colsize());
		gsl_matrix* result_err = gsl_matrix_alloc(TwoPointFunction.rowsize(), TwoPointFunction.colsize());
		gsl_matrix* result_err_rel = gsl_matrix_alloc(TwoPointFunction.rowsize(), TwoPointFunction.colsize());

		for (size_type i = 0; i < result->size1; ++i)
			{
				for (size_type j = 0; j < result->size2; ++j)
				{
					number_type real = TwoPointFunction.get_real(i, j);
					number_type imag = TwoPointFunction.get_imag(i, j);
					//gsl_matrix_set(moduli, i, j, sqrt(real*real+imag*imag));
					gsl_matrix_set(result, i, j, real);
					gsl_matrix_set(result_err, i, j, TwoPointFunction_err.get_real(i, j));
					gsl_matrix_set(result_err_rel, i, j, TwoPointFunction_err.get_real(i, j)/sqrt(real*real+imag*imag));
				}
			}
			std::cout << "Maximal absolute error: " << gsl_matrix_max(result_err) << "\n"; 
			number_type maximal_modulus = std::max(gsl_matrix_max(result), -gsl_matrix_min(result));
			number_type current_max_rel_error = gsl_matrix_max(result_err)/maximal_modulus;
			std::cout << "In relation to maximal value: " << current_max_rel_error << "\n"; 
			if (current_max_rel_error > max_rel_error)
			{
				max_rel_error = current_max_rel_error;
			}
			to_file(outfile_background, result);
			to_file(outfile_background_err, result_err);
			gsl_matrix_free(result);
			gsl_matrix_free(result_err);
			gsl_matrix_free(result_err_rel);
	}

	std::cout << "Maximal relative error found: " << max_rel_error << "\n";

	max_rel_error = -1.0;

	for (int m = 0; m < 5; ++m)
	{
		for (int c = 1; c < classes.size(); ++c)
		{
			std::cout << "m: " << m << " class: " << classes[c] << "%" << "\n";
			// Create outfile name
			std::string centrality_class = std::to_string(int(classes[c-1])) + "-" + std::to_string(int(classes[c]));
			std::string outfile_modulus = "output/two_point_random_" + centrality_class;
			outfile_modulus += "_m";
			outfile_modulus += std::to_string(m);
			outfile_modulus += "_real";
			std::string outfile_modulus_err = outfile_modulus;
			outfile_modulus_err += "_error.txt";
			outfile_modulus += ".txt";
			size_type lMax = 10;

			complex_matrix<number_type> TwoPointFunction(lMax, lMax);
			complex_matrix<number_type> TwoPointFunction_err(lMax, lMax);
			PbPb.getTwoPointFunction_fixed_m(m, c, TwoPointFunction, TwoPointFunction_err, start);


			// save moduli of coeffs in text file 
			gsl_matrix* moduli = gsl_matrix_alloc(TwoPointFunction.rowsize(), TwoPointFunction.colsize());
			gsl_matrix* moduli_error = gsl_matrix_alloc(TwoPointFunction.rowsize(), TwoPointFunction.colsize());
			gsl_matrix* moduli_error_rel = gsl_matrix_alloc(TwoPointFunction.rowsize(), TwoPointFunction.colsize());
			for (size_type i = 0; i < moduli->size1; ++i)
			{
				for (size_type j = 0; j < moduli->size2; ++j)
				{
					number_type real = TwoPointFunction.get_real(i, j);
					number_type imag = TwoPointFunction.get_imag(i, j);
					//gsl_matrix_set(moduli, i, j, sqrt(real*real+imag*imag));
					gsl_matrix_set(moduli, i, j, real);
					gsl_matrix_set(moduli_error, i, j, TwoPointFunction_err.get_real(i, j));
					gsl_matrix_set(moduli_error_rel, i, j, TwoPointFunction_err.get_real(i, j)/sqrt(real*real+imag*imag));
				}
			}
			std::cout << "Maximal absolute error: " << gsl_matrix_max(moduli_error) << "\n"; 
			number_type maximal_modulus = std::max(gsl_matrix_max(moduli), -gsl_matrix_min(moduli));
			number_type current_max_rel_error = gsl_matrix_max(moduli_error)/maximal_modulus;
			std::cout << "In relation to maximal value: " << current_max_rel_error << "\n"; 
			if (current_max_rel_error > max_rel_error)
			{
				max_rel_error = current_max_rel_error;
			}
			to_file(outfile_modulus, moduli);
			to_file(outfile_modulus_err, moduli_error);
			gsl_matrix_free(moduli);
			gsl_matrix_free(moduli_error);
			gsl_matrix_free(moduli_error_rel);

			// save phase of coeffs in text file
			std::string outfile_phase = "output/two_point_random_" + centrality_class;
			outfile_phase += "_m";
			outfile_phase += std::to_string(m);
			outfile_phase += "_phase";
			outfile_phase += ".txt";
			gsl_matrix* phases = gsl_matrix_alloc(TwoPointFunction.rowsize(), TwoPointFunction.colsize());
			for (size_type i = 0; i < phases->size1; ++i)
			{
				for (size_type j = 0; j < phases->size2; ++j)
				{
					number_type real = TwoPointFunction.get_real(i, j);
					number_type imag = TwoPointFunction.get_imag(i, j);
					gsl_matrix_set(phases, i, j, atan2(imag, real));
				}
			} 
			to_file(outfile_phase, phases);
			gsl_matrix_free(phases);

			

		}
	}
	


	current_time = std::time(nullptr);
	std::cout << current_time-start << "s: " << "Output data saved.\n"; 

	std::cout << "Maximal relative error found: " << max_rel_error << "\n";




	return 0;
}