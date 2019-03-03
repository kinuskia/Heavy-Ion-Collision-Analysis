
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

	//PbPb.normalize(1); // normalize events so that integral = 1


	PbPb.centralize(); // shift data so that barycentre at origin

	PbPb.getReactionPlane("output/angles.txt"); // compute reaction plane angles




	// define the respective multiplicity limits for specific centrality classes
	std::vector<number_type> classes(6);
	classes[0] = 0;
	classes[1] = 5;
	classes[2] = 10;
	classes[3] = 20; // this generates the classes 0-5%, 5-10, 10-20, 20-30, ..., 80-90, 90-100
	classes[4] = 30;
	classes[5] = 40;
	PbPb.get_percentiles(classes);

	std::time_t current_time = std::time(nullptr);
	std::cout << current_time-start << "s: " << "Data has been read in and normalized.\n"; 


	// print ensemble averaged profiles
	PbPb.print_averaged_profiles("output/profiles_averaged");

	// set interpolation method: gsl_interp2d_bicubic or gsl_interp2d_bilinear
	const gsl_interp2d_type* xy_interpolation_method = gsl_interp2d_bicubic;
	PbPb.initialize_xy_interpolation(xy_interpolation_method); // initialize objects needed for interpolation

	

	/*
	compute expectation values <e_0l> for each centrality class (m=0 for non-vanishing values)
	*/

	const gsl_interp_type* r_interpolation_method = gsl_interp_cspline;
	//PbPb.initialize_n_point_evaluations(r_interpolation_method, start);
	/* Fourier-decompose the profiles */
	PbPb.decompose_azimuthal();
	PbPb.getNormalizations();

	// print <e_m(r)> for all centrality classes
	for (size_type m = 0; m < 5; ++m)
	{
		for (size_type c = 1; c < classes.size(); ++c)
		{
			std::string filename_mean = "output/e_m_mean_m";
			filename_mean += std::to_string(m);
			filename_mean += "_";
			std::string centrality_class = std::to_string(int(classes[c-1])) + "-" + std::to_string(int(classes[c]));
			filename_mean += centrality_class;
			filename_mean += ".txt";
			std::cout << "m: " << m << " c: " << c << "\n";
			PbPb.print_e_m(filename_mean, m, c);
		}
	}
	

	PbPb.initialize_r_interpolation(r_interpolation_method);

	// print weighting functions W
	PbPb.print_W("output/weight_functions.txt", 200);



	// print maps rho
	PbPb.print_rho("output/rhos.txt", 200);

	//std::cout << PbPb.integ_test() << "\n";

	// compute expectation value of background coefficient
	std::vector<number_type> eps_00(0);
	std::vector<number_type> eps_00_err(0);
	for (int c = 1; c < classes.size(); ++c)
	{
		std::cout << " class: " << classes[c] << "%" << "\n";
		number_type coeff_mean;
		number_type coeff_err;
		PbPb.getOnePointFunction_00(c, coeff_mean, coeff_err, start);
		eps_00.push_back(coeff_mean);
		eps_00_err.push_back(coeff_err);
	}
	std::vector<std::vector<number_type>> eps_00_data(2);
	eps_00_data[0] = eps_00;
	eps_00_data[1] = eps_00_err;
	to_file("output/background_coeffs.txt", eps_00_data);


	for (int c = 1; c < classes.size(); ++c)
	{
		std::cout << " class: " << classes[c] << "%" << "\n";
		// Create outfile names
		std::string centrality_class = std::to_string(int(classes[c-1])) + "-" + std::to_string(int(classes[c]));
		std::string outfile_modulus = "output/one_point_" + centrality_class;
		//outfile_modulus += "_modulus";
		outfile_modulus += ".txt";
		size_type lMax = 10;
		size_type mMax = 5;

		complex_matrix<number_type> OnePointFunction(mMax+1, lMax);
		complex_matrix<number_type> OnePointFunction_err(mMax+1, lMax);

		PbPb.getOnePointFunction(c, OnePointFunction, OnePointFunction_err, start);
		
		// save real parts of coeffs in text file
		gsl_matrix* result = gsl_matrix_alloc(mMax+1, lMax); 
		for (size_type i = 0; i < result->size1; ++i)
		{
			for (size_type j = 0; j < result->size2; ++j)
			{
				number_type real = OnePointFunction.get_real(i, j);
				//number_type imag = OnePointFunction.get_imag(i, j);
				gsl_matrix_set(result, i, j, real);	
			}
		}
 
		to_file(outfile_modulus, result);

		// save error of coeffs in text file
		std::string outfile_phase = "output/one_point_" + centrality_class;
		outfile_phase += "_error";
		outfile_phase += ".txt";
		gsl_matrix* result_error = gsl_matrix_alloc(mMax+1, lMax);
		for (size_type i = 0; i < result_error->size1; ++i)
			for(size_type j = 0; j < result_error->size2; ++j)
			{
				number_type real = OnePointFunction_err.get_real(i, j);
				//number_type imag = OnePointFunction.get_imag(i, j);
				gsl_matrix_set(result_error, i, j, real);	
			}


		to_file(outfile_phase, result_error); 

	}

	// Compute clm values
	size_type lMax = 20; 
	size_type mMax = 10;
	PbPb.get_Bessel_deriv_zeros(mMax+1, lMax+1);
	gsl_matrix* clm = gsl_matrix_alloc(lMax, mMax +1);
	for (size_type l = 1; l<= lMax ; ++l)
	{
		for (size_type m = 0; m <= mMax; ++m)
		{
			gsl_matrix_set(clm, l-1, m, PbPb.c(m, l));
		}
	}
	to_file("output/clm.txt", clm);

	PbPb.print_Bessel_deriv_zeros(10, 5, "output/bessel_d_0.txt");






	current_time = std::time(nullptr);
	std::cout << current_time-start << "s: " << "Output data saved.\n"; 




	return 0;
}