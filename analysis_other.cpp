
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
	std::vector<number_type> classes(22);
	classes[0] = 0;
	classes[1] = 1;
	classes[2] = 2;
	classes[3] = 3; 
	classes[4] = 4;
	classes[5] = 5;
	classes[6] = 6;
	classes[7] = 7;
	classes[8] = 8;
	classes[9] = 9; 
	classes[10] = 10;
	classes[11] = 11;
	classes[12] = 12;
	classes[13] = 13;
	classes[14] = 14; 
	classes[15] = 15;
	classes[16] = 16;
	classes[17] = 17;
	classes[18] = 18;
	classes[19] = 19;
	classes[20] = 20; 
	classes[21] = 21;
	PbPb.get_percentiles(classes);

	std::time_t current_time = std::time(nullptr);
	std::cout << current_time-start << "s: " << "Data has been read in. \n"; 


	// print ensemble averaged profiles
	//PbPb.print_averaged_profiles("output/profiles_averaged");

	// set interpolation method: gsl_interp2d_bicubic or gsl_interp2d_bilinear
	const gsl_interp2d_type* xy_interpolation_method = gsl_interp2d_bicubic;
	PbPb.initialize_xy_interpolation(xy_interpolation_method); // initialize objects needed for interpolation

	

	/*
	compute expectation values <e_0l> for each centrality class (m=0 for non-vanishing values)
	*/

	
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
	
	const gsl_interp_type* r_interpolation_method = gsl_interp_cspline;
	PbPb.initialize_r_interpolation(r_interpolation_method);

	// print weighting functions W
	PbPb.print_W("output/weight_functions.txt", 200);



	// print maps rho
	PbPb.print_rho("output/rhos.txt", 200);

	//std::cout << PbPb.integ_test() << "\n";


	// Print file with the following columns: impact parameter, number of participants, multiplicity
	PbPb.collision_specs_to_file("output/collision_specs.txt");
	
	// Print collision specs file for each centrality class
	
	for (size_type c = 1; c < classes.size(); ++c)
	{
		std::string filename = "output/collision_specs_";
		std::string centrality_class = std::to_string(int(classes[c-1])) + "-" + std::to_string(int(classes[c]));
		filename += centrality_class;
		filename += ".txt";

		PbPb.collision_specs_to_file(filename, c);
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