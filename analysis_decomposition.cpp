
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
	PbPb.get_Bessel_deriv_zeros("auxiliary/bessel_zeros.txt", 11, 5);

	PbPb.read_in(filename, fileformat, n_files); // read in Trento event files

	PbPb.normalize(1); // normalize events so that integral = 1


	PbPb.centralize(); // shift data so that barycentre at origin

	// set interpolation method: gsl_interp2d_bicubic or gsl_interp2d_bilinear
	const gsl_interp2d_type* xy_interpolation_method = gsl_interp2d_bicubic;
	PbPb.initialize_xy_interpolation(xy_interpolation_method); // initialize objects needed for interpolation


	/*
	goal generate text file with following columns
	m 	l 	e_ml	err_ml
	*/
	// Create outfile name
	std::string outfile = "output/decomposition";
	outfile += impact_parameter;
	outfile += ".txt";
	size_type mMax = 5;
	size_type lMax = 20;
	complex_matrix<number_type> Fourier_Bessel_coeffs_mean(mMax+1, lMax);
	complex_matrix<number_type> Fourier_Bessel_coeffs_err(mMax+1, lMax);
	const gsl_interp_type* r_interpolation_method = gsl_interp_cspline;
	PbPb.FourierBesselDecompose(Fourier_Bessel_coeffs_mean, Fourier_Bessel_coeffs_err, r_interpolation_method);

	Fourier_Bessel_coeffs_mean.print();

	Fourier_Bessel_coeffs_err.print();

	// save moduli of coeffs in text file 
	gsl_matrix* moduli = gsl_matrix_alloc(Fourier_Bessel_coeffs_mean.rowsize(), Fourier_Bessel_coeffs_mean.colsize());
	for (size_type i = 0; i < moduli->size1; ++i)
	{
		for (size_type j = 0; j < moduli->size2; ++j)
		{
			number_type real = Fourier_Bessel_coeffs_mean.get_real(i, j);
			number_type imag = Fourier_Bessel_coeffs_mean.get_imag(i, j);
			gsl_matrix_set(moduli, i, j, sqrt(real*real+imag*imag));
		}
	} 
	to_file(outfile, moduli);
	gsl_matrix_free(moduli);

	PbPb.free();




	return 0;
}