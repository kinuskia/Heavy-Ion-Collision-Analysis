#ifndef FBDECOMPOSITION_HPP
#define FBDECOMPOSITION_HPP

#include <iostream>
#include <cmath>

#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_matrix.h>

#include <gsl/gsl_fft_real.h>
#include <algorithm>

#include <gsl/gsl_sf_bessel.h>
#include "../auxiliary/bessel_deriv_zero.hpp"

//#include "../complex_matrix.hpp"


template<class REAL>
class FBDecomposition
{
public:
	typedef std::size_t size_type;
	typedef REAL number_type;

	//constructor
	FBDecomposition(
	  Model<number_type> model
	, number_type rMax
	 )
	: model_(model)
	, rMax_(rMax)
	, pi_(3.14159265359)
	, normalization_(1.)
	, Nm_(128)
	, Nr_(50)
	, N_discret_(25)
	{

	}

	// method for initializing the object
	void initialize()
	{
		// Compute one-point normalization
		DetermineNormalization();

		// Initialize W(r)
		initialize_W();

		// Initialize rho(r)
		initialize_rho();
	}

	// getter for rMax_
	number_type get_rMax() const
	{
		return rMax_;
	}

	// setter for N_discret_
	void set_N_discret(size_type N)
	{
		N_discret_ = N;
	}

	// integrate one-point function over phi at fixed r
	// (only declaration, it is defined outside this class, see below)
	number_type integ_one_phi(number_type r, number_type phi_lower, number_type phi_upper);

	// integrate (phi-integrated) one-point function over r
	// so for total integral over one-point function, first apply "integ_one_phi", then this function.
	// (only declaration, it is defined outside this class, see below)
	number_type integ_one_r(number_type r_lower, number_type r_upper);


	// Getter for one-Point correlation function in polar coordinates
	number_type OnePoint(number_type r, number_type phi)
	{
		number_type x = r*cos(phi);
		number_type y = r*sin(phi);

		return model_.OnePoint(x, y)/normalization_;
	}

	// Getter for two-Point correlation function in polar coordinates
	number_type TwoPoint(number_type r1, number_type phi1, number_type r2, number_type phi2)
	{
		number_type x1 = r1*cos(phi1);
		number_type y1 = r1*sin(phi1);
		number_type x2 = r2*cos(phi2);
		number_type y2 = r2*sin(phi2);

		return model_.TwoPoint(x1, y1, x2, y2);
	}

	// compute total integral of one point function = normalization
	void DetermineNormalization()
	{
		normalization_ = integ_one_r(0, rMax_);
	}

	// initialize W(r) function
	void initialize_W()
	{
		// Compute certain grid points of W and generate a spline
		number_type* W_sites;
		number_type* r_sites;
		W_sites = new number_type[Nr_];
		r_sites = new number_type[Nr_];
		for (size_type i = 0; i < Nr_; ++i)
		{
			number_type r = rMax_*i/Nr_;
			r_sites[i] = r;
			W_sites[i] = integ_one_phi(r, 0, 2.*pi_)/2;
		}
		const gsl_interp_type* interpolator = gsl_interp_cspline;
		gsl_spline* spline = gsl_spline_alloc(interpolator, Nr_);
		gsl_interp_accel* acc = gsl_interp_accel_alloc();

		W_interpolator_ = interpolator;
		W_spline_ = spline;
		W_acc_ = acc;
		gsl_spline_init(W_spline_, r_sites, W_sites, Nr_);
	}

	// getter for W(r) function
	number_type W(number_type r)
	{
		return gsl_spline_eval(W_spline_, r, W_acc_);
	}

	// integrate W(r)r over r
	// definition outside class
	number_type integ_W_r(number_type r_lower, number_type r_upper);

	// initialize rho(r) function
	void initialize_rho()
	{
		// Compute certain grid points of W and generate a spline
		number_type* rho_sites;
		number_type* r_sites;
		rho_sites = new number_type[Nr_];
		r_sites = new number_type[Nr_];
		for (size_type i = 0; i < Nr_; ++i)
		{
			number_type r = rMax_*i/Nr_;
			r_sites[i] = r;
			rho_sites[i] = sqrt(integ_W_r(0, r)*2);
		}
		const gsl_interp_type* interpolator = gsl_interp_cspline;
		gsl_spline* spline = gsl_spline_alloc(interpolator, Nr_);
		gsl_interp_accel* acc = gsl_interp_accel_alloc();

		rho_interpolator_ = interpolator;
		rho_spline_ = spline;
		rho_acc_ = acc;
		gsl_spline_init(rho_spline_, r_sites, rho_sites, Nr_);
	}

	// getter for rho(r) function
	number_type rho(number_type r)
	{
		return gsl_spline_eval(rho_spline_, r, rho_acc_);
	}

	// Compute zeros of first derivative of the Bessel function J (l'th zero crossing of J'_m)
	void get_Bessel_deriv_zeros(size_type m, size_type l)
	{
		bessel_deriv_zeros_ = gsl_matrix_alloc(m, l);
		number_type eps_rel = 1.e-15;
		size_type max_iter = 100;
		for (size_type i = 0; i < bessel_deriv_zeros_->size1; ++i)
		{
			for (size_type j = 0; j < bessel_deriv_zeros_->size2; ++j)
			{
				gsl_matrix_set(bessel_deriv_zeros_, i, j, find_bessel_deriv_root(i, j+1, eps_rel, max_iter));
			}
		}
	}

	number_type Bessel_zero(int m, size_type l)
	{
		if (m < 0) // take absolute value of m
		{
			m = -m;
		}

		return gsl_matrix_get(bessel_deriv_zeros_, m, l-1);
		//return gsl_sf_bessel_zero_Jnu(m, l);
	}

	// Evaluate appropriate Bessel function
	number_type Bessel(int m, size_type l, number_type r)
	{
		return gsl_sf_bessel_Jn(m, Bessel_zero(m, l)*rho(r));
	}

	// Evaluate correct weighting function in the r integration
	number_type weight(int m, size_type l, number_type r)
	{
		number_type result;
		if (m == 0 && l == 1)
		{
			result =  r; 
		}
		else if (m == 0 && l > 1)
		{
			result = r*Bessel(m, l-1, r);
		}
		else
		{
			result = r*Bessel(m, l, r);
		}

		return result;
	}

	// evaluate Bessel coefficient c_ml (if one uses J' zeros)
	number_type c (int m, size_type l)
	{
		number_type zero = Bessel_zero(m, l);
		number_type result = gsl_sf_bessel_Jn(m, zero)*gsl_sf_bessel_Jn(m, zero)*(zero*zero - m*m)/2./zero/zero;
		return result; 
	}

	// integrate two-point function over r2 at fixed phi1, phi2, r1, m2, l2
	// and Besselfunction in r2
	number_type integ_r2(int m2, int l2, number_type phi1, number_type phi2, number_type r1);

	// integrate two-point function over r1 at fixed phi1, phi2, m1, l1
	// and Besselfunction in r1
	number_type integ_r1(int m1, int l1, int m2, int l2, number_type phi1, number_type phi2);

	
	// my own integration routine for the r1, r2 integration
	number_type integ_r(int m1, int l1, int m2, int l2, number_type phi1, number_type phi2)
	{
		number_type result = 0;
		size_type N = N_discret_;
		number_type width = (rMax_-0.2)/N;

		// compute integral over r1
		// Use trapezoidal rule
		for (size_type i = 0; i <= N; ++i) // loop over r1
		{
			number_type r1 = (rMax_-0.2)*i/N;
			//std::cout << "r1: " << r1 << "\n";

			// compute integral over r2 at fixed r1
			number_type integral_r2 = 0;
			for (size_type j = 0; j < N; ++j) // one fewer grid point and offset so that f is never evaluated at the same radius as r1
			{
				number_type r2 = (rMax_-0.2)*(0.5+j)/N;
				number_type f = weight(m1, l1, r1)*weight(m2, l2, r2)*TwoPoint(r1, phi1, r2, phi2);
				//std::cout << "phi1: " << phi1 << " phi2: " << phi2  <<  " r2: " << r2 << " f: " << f << "\n";
				// Use trapezoidal rule
				if ((j == 0) || (j == (N-1)))
				{
					integral_r2 += f/2;
				}
				else
				{
					integral_r2 += f;
				}
			}
			integral_r2 *= width;

			// use integral_r2 result as input for trapezoidal rule in r1-direction
			if ((i == 0) || (i == N))
			{
				result += integral_r2/2;
			}
			else
			{
				result += integral_r2;
			}
		}
		result *= width;

		//std::cout << result << "\n";


		return result;
	}

	// compute two-mode correlator <e_l1^(m1) e_l2^(m2)>
	number_type TwoMode(int m1, int l1, int m2, int l2)
	{
		get_Bessel_deriv_zeros(std::max(abs(m1), abs(m2))+1, std::max(l1, l2));
		//Compute FFT of a function E_l1l2^(m1m2)(phi1, phi2)
		
		// Step 1: Compute FFT with respect to phi1 at fixed phi2
		// For each phi2, the FFT with m=m1 is saved in a vector
		number_type* FFT_real;
		number_type* FFT_imag;
		FFT_real = new number_type[Nm_];
		FFT_imag = new number_type[Nm_];
		
		// Loop over phi2
		for (size_type j = 0; j < Nm_; ++j)
		{
			// Do FFT with respect to phi1 at fixed phi2
			number_type phi2 = 2.* pi_ * j / Nm_;
			number_type* fft;
			fft = new number_type[Nm_];

			for (size_type i = 0; i < Nm_; ++i)
			{
				number_type phi1 = 2. * pi_ * i / Nm_;
			
				fft[i] = integ_r(m1, l1, m2, l2, phi1, phi2);
				// 	if (i == 0)
				// {
				// 	std::cout << "m1: " << m1 << " ";
				// 	std::cout << "l1: " << l1 << " ";
				// 	std::cout << "m2: " << m2 << " ";
				// 	std::cout << "l2: " << l2 << " ";
				// 	std::cout << fft[i] << "\n";
				// }
			}

			gsl_fft_real_radix2_transform(fft, 1, Nm_);

			// save result
			if (m1 >= 0)
			{
				FFT_real[j] = fft[m1]/Nm_;
				if (m1 == 0)
				{
					FFT_imag[j] = 0;
				}
				else
				{
					FFT_imag[j] = fft[Nm_-m1]/Nm_;
				}
			}
			else
			{
				FFT_real[j] = fft[-m1]/Nm_;
				FFT_imag[j] = -fft[Nm_+m1]/Nm_;
			}
		}

		// Step 2: Compute FFT with respect to phi2
		gsl_fft_real_radix2_transform(FFT_real, 1, Nm_);
		gsl_fft_real_radix2_transform(FFT_imag, 1, Nm_);
		number_type result;
		if (m2 == 0)
		{
			result = (FFT_real[m2])/Nm_;
		}
		else if (m2 > 0)
		{
			result = (FFT_real[m2] - FFT_imag[Nm_-m2])/Nm_;
		}
		else
		{
			result = (FFT_real[-m2] + FFT_imag[Nm_+m2])/Nm_;
		}



		// Step 3: Scale result with the right c_l^(m)
		if (m1 == 0 && l1 == 1)
		{
			result /= 0.5;
		}
		else if (m1 == 0 && l1 > 1)
		{
			result /= c(m1, l1-1);
		}
		else
		{
			result /= c(m1, l1);
		}

		if (m2 == 0 && l2 == 1)
		{
			result /= 0.5;
		}
		else if (m2 == 0 && l2 > 1)
		{
			result /= c(m2, l2-1);
		}
		else
		{
			result /= c(m2, l2);
		}

		return result;
	}

	// compute two-mode correlator <e_l1^(m1) e_l2^(m2)>
	number_type TwoMode_imag(int m1, int l1, int m2, int l2)
	{
		get_Bessel_deriv_zeros(std::max(abs(m1), abs(m2))+1, std::max(l1, l2));
		//Compute FFT of a function E_l1l2^(m1m2)(phi1, phi2)
		
		// Step 1: Compute FFT with respect to phi1 at fixed phi2
		// For each phi2, the FFT with m=m1 is saved in a vector
		number_type* FFT_real;
		number_type* FFT_imag;
		FFT_real = new number_type[Nm_];
		FFT_imag = new number_type[Nm_];
		
		// Loop over phi2
		for (size_type j = 0; j < Nm_; ++j)
		{
			// Do FFT with respect to phi1 at fixed phi2
			number_type phi2 = 2.* pi_ * j / Nm_;
			number_type* fft;
			fft = new number_type[Nm_];

			for (size_type i = 0; i < Nm_; ++i)
			{
				number_type phi1 = 2. * pi_ * i / Nm_;
				fft[i] = integ_r1(m1, l1, m2, l2, phi1, phi2);
			}

			gsl_fft_real_radix2_transform(fft, 1, Nm_);

			// save result 
			if (m1 >= 0)
			{
				FFT_real[j] = fft[m1]/Nm_;
				if (m1 == 0)
				{
					FFT_imag[j] = 0;
				}
				else
				{
					FFT_imag[j] = fft[Nm_-m1]/Nm_;
				}
			}
			else
			{
				FFT_real[j] = fft[-m1]/Nm_;
				FFT_imag[j] = -fft[Nm_+m1]/Nm_;
			}
		}

		// Step 2: Compute FFT with respect to phi2
		gsl_fft_real_radix2_transform(FFT_real, 1, Nm_);
		gsl_fft_real_radix2_transform(FFT_imag, 1, Nm_);
		number_type result;
		if (m2 == 0)
		{
			result = (FFT_imag[m2])/Nm_;
		}
		else if (m2 > 0)
		{
			result = (FFT_real[Nm_-m2] + FFT_imag[m2])/Nm_;
		}
		else
		{
			result = (-FFT_real[Nm_+m2] + FFT_imag[-m2])/Nm_;
		}



		// Step 3: Scale result with the right c_l^(m)
		if (m1 == 0 && l1 == 1)
		{
			result /= 0.5;
		}
		else if (m1 == 0 && l1 > 1)
		{
			result /= c(m1, l1-1);
		}
		else
		{
			result /= c(m1, l1);
		}

		if (m2 == 0 && l2 == 1)
		{
			result /= 0.5;
		}
		else if (m2 == 0 && l2 > 1)
		{
			result /= c(m2, l2-1);
		}
		else
		{
			result /= c(m2, l2);
		}

		return result;
	}

	// 2D Fourier test
	number_type Fourier2D(int m1, int m2)
	{
		
		// Step 1: Compute FFT with respect to phi1 at fixed phi2
		// For each phi2, the FFT with m=m1 is saved in a vector
		number_type* FFT_real;
		number_type* FFT_imag;
		FFT_real = new number_type[Nm_];
		FFT_imag = new number_type[Nm_];
		
		// Loop over phi2
		for (size_type j = 0; j < Nm_; ++j)
		{
			// Do FFT with respect to phi1 at fixed phi2
			number_type phi2 = 2.* pi_ * j / Nm_;
			number_type* fft;
			fft = new number_type[Nm_];

			for (size_type i = 0; i < Nm_; ++i)
			{
				number_type phi1 = 2. * pi_ * i / Nm_;
				fft[i] = (phi1-phi2)*(phi1-phi2);
			}

			gsl_fft_real_radix2_transform(fft, 1, Nm_);

			// save result (only real part)
			if (m1 >= 0)
			{
				FFT_real[j] = fft[m1]/Nm_;
				if (m1 == 0)
				{
					FFT_imag[j] = 0;
				}
				else
				{
					FFT_imag[j] = fft[Nm_-m1]/Nm_;
				}
			}
			else
			{
				FFT_real[j] = fft[-m1]/Nm_;
				FFT_imag[j] = -fft[Nm_+m1]/Nm_;
			}
		}

		// Step 2: Compute FFT with respect to phi2
		gsl_fft_real_radix2_transform(FFT_real, 1, Nm_);
		gsl_fft_real_radix2_transform(FFT_imag, 1, Nm_);
		number_type result;
		if (m2 == 0)
		{
			result = (FFT_real[m2])/Nm_;
		}
		else if (m2 > 0)
		{
			result = (FFT_real[m2] - FFT_imag[Nm_-m2])/Nm_;
		}
		else
		{
			result = (FFT_real[-m2] + FFT_imag[Nm_+m2])/Nm_;
		}

		return result;
	}

	number_type pi_;

private:
	Model<number_type> model_;
	number_type rMax_;
	size_type Nr_;
	size_type Nm_;
	number_type normalization_;

	size_type N_discret_;

	const gsl_interp_type* W_interpolator_;
	gsl_spline* W_spline_;
	gsl_interp_accel* W_acc_;

	const gsl_interp_type* rho_interpolator_;
	gsl_spline* rho_spline_;
	gsl_interp_accel* rho_acc_;

	gsl_matrix* bessel_deriv_zeros_;
};

















// struct to feed radius to phi integration function
template<typename number_type>
struct integ_one_phi_params
{
	number_type radius;
	FBDecomposition<number_type>* pt_FBDecomposition;
};

// gsl function to call to integrate over one-point function with respect to phi at fixed r
double integ_one_phi_wrapper(double phi, void* params)
{
	integ_one_phi_params<double>* params_ = (integ_one_phi_params<double>*) params;
	double r = params_->radius;
	return params_->pt_FBDecomposition->OnePoint(r, phi);
}

// integrate one-point function over phi at fixed r
template<class number_type>
number_type FBDecomposition<number_type>::integ_one_phi(number_type r, number_type phi_lower, number_type phi_upper)
{
	//typedef std::size_t size_type;
	gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);
	number_type integral;
	number_type error;
	integ_one_phi_params<number_type> params = {r, this};

	gsl_function F;
	F.function = &integ_one_phi_wrapper;
	F.params = &params;
	gsl_integration_qags(&F, phi_lower, phi_upper, 1e-5, 1e-5, 1000, w, &integral, &error);
	//std::cout << " Radius: " << r << " Index: " << k <<  " Integral: " << integral << " Error abs: " << error << " Error rel: " << error/integral << "\n"; 
	gsl_integration_workspace_free(w);
	return integral;
}

// struct to feed class to r integration function
template<typename number_type>
struct integ_one_r_params
{
	FBDecomposition<number_type>* pt_FBDecomposition;
};

// gsl function to call to integrate over one-point function with respect to phi at fixed r
double integ_one_r_wrapper(double r, void* params)
{
	integ_one_r_params<double>* params_ = (integ_one_r_params<double>*) params;
	double pi = params_->pt_FBDecomposition->pi_;
	return (params_->pt_FBDecomposition->integ_one_phi(r, 0, 2.*pi))*r; // times r because of Jacobian
}

// integrate (phi-integrated) one-point function over r
// so for total integral over one-point function, first apply "integ_one_phi", then this function.
template<class number_type>
number_type FBDecomposition<number_type>::integ_one_r(number_type r_lower, number_type r_upper)
{
	//typedef std::size_t size_type;
	gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);
	number_type integral;
	number_type error;
	integ_one_r_params<number_type> params = {this};

	gsl_function F;
	F.function = &integ_one_r_wrapper;
	F.params = &params;
	gsl_integration_qags(&F, r_lower, r_upper, 1e-5, 1e-5, 1000, w, &integral, &error);
	//std::cout << " Radius: " << r << " Index: " << k <<  " Integral: " << integral << " Error abs: " << error << " Error rel: " << error/integral << "\n"; 
	gsl_integration_workspace_free(w);
	return integral;
}

// struct for integ_W_r
template<typename number_type>
struct integ_W_r_params
{
	FBDecomposition<number_type>* pt_FBDecomposition;
};

// gsl function to call to integrate over W(r)r
double integ_W_r_wrapper(double r, void* params)
{
	integ_W_r_params<double>* params_ = (integ_W_r_params<double>*) params;

	return (params_->pt_FBDecomposition->W(r))*r;
}

// integrate W(r)r over r
template<class number_type>
number_type FBDecomposition<number_type>::integ_W_r(number_type r_lower, number_type r_upper)
{
	//typedef std::size_t size_type;
	gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);
	number_type integral;
	number_type error;
	integ_W_r_params<number_type> params = {this};

	gsl_function F;
	F.function = &integ_W_r_wrapper;
	F.params = &params;
	gsl_integration_qags(&F, r_lower, r_upper, 1e-5, 1e-5, 1000, w, &integral, &error);
	//std::cout << " Radius: " << r << " Index: " << k <<  " Integral: " << integral << " Error abs: " << error << " Error rel: " << error/integral << "\n"; 
	gsl_integration_workspace_free(w);
	return integral;
}

// struct for integ_r2
template<typename size_type, typename number_type>
struct integ_r2_params
{
	size_type m2;
	size_type l2;
	number_type phi1;
	number_type phi2;
	number_type r1;
	FBDecomposition<number_type>* pt_FBDecomposition;
};

// gsl function to call to integrate over r2
double integ_r2_wrapper(double r, void* params)
{
	typedef int size_type;
	typedef double number_type;
	integ_r2_params<size_type, number_type>* params_ = (integ_r2_params<size_type, number_type>*) params;
	size_type m2_ = params_->m2;
	size_type l2_ = params_->l2;
	number_type phi1_ = params_->phi1;
	number_type phi2_ = params_->phi2;
	number_type r1_ = params_->r1;
	return (params_->pt_FBDecomposition->weight(m2_, l2_, r))*(params_->pt_FBDecomposition->TwoPoint(r1_, phi1_, r, phi2_));
}


// integrate two-point function over r2 at fixed phi1, phi2, r1, m2, l2
// and Besselfunction in r2
template<class number_type>
number_type FBDecomposition<number_type>::integ_r2(int m2, int l2, number_type phi1, number_type phi2, number_type r1)
{
	typedef int size_type;
	gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);
	number_type integral;
	number_type error;
	integ_r2_params<size_type, number_type> params = {m2, l2, phi1, phi2, r1, this};

	gsl_function F;
	F.function = &integ_r2_wrapper;
	F.params = &params;
	gsl_integration_qags(&F, 0, (this->get_rMax()-0.2), 1e-5, 1e-5, 1000, w, &integral, &error);
	//std::cout << " Radius: " << r << " Index: " << k <<  " Integral: " << integral << " Error abs: " << error << " Error rel: " << error/integral << "\n"; 
	gsl_integration_workspace_free(w);
	return integral;
}

// struct for integ_r1
template<typename size_type, typename number_type>
struct integ_r1_params
{
	size_type m1;
	size_type l1;
	size_type m2;
	size_type l2;
	number_type phi1;
	number_type phi2;
	FBDecomposition<number_type>* pt_FBDecomposition;
};

// gsl function to call to integrate over r1
double integ_r1_wrapper(double r, void* params)
{
	typedef int size_type;
	typedef double number_type;
	integ_r1_params<size_type, number_type>* params_ = (integ_r1_params<size_type, number_type>*) params;
	size_type m1_ = params_->m1;
	size_type l1_ = params_->l1;
	size_type m2_ = params_->m2;
	size_type l2_ = params_->l2;
	number_type phi1_ = params_->phi1;
	number_type phi2_ = params_->phi2;
	return (params_->pt_FBDecomposition->weight(m1_, l1_, r))*(params_->pt_FBDecomposition->integ_r2(m2_, l2_, phi1_, phi2_, r));
}


// integrate two-point function over r1 at fixed phi1, phi2, m2, l2
// and Besselfunction in r1
template<class number_type>
number_type FBDecomposition<number_type>::integ_r1(int m1, int l1, int m2, int l2, number_type phi1, number_type phi2)
{
	typedef int size_type;
	gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);
	number_type integral;
	number_type error;
	integ_r1_params<size_type, number_type> params = {m1, l1, m2, l2, phi1, phi2, this};

	gsl_function F;
	F.function = &integ_r1_wrapper;
	F.params = &params;
	gsl_integration_qags(&F, 0, (this->get_rMax()-0.2), 1e-5, 1e-5, 1000, w, &integral, &error);
	//std::cout << " Radius: " << r << " Index: " << k <<  " Integral: " << integral << " Error abs: " << error << " Error rel: " << error/integral << "\n"; 
	gsl_integration_workspace_free(w);
	return integral;
}

#endif