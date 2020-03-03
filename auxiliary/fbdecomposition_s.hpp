#ifndef FBDECOMPOSITION_S_HPP
#define FBDECOMPOSITION_S_HPP

#include <iostream>
#include <cmath>

#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_matrix.h>

#include <gsl/gsl_fft_real.h>
#include <algorithm>

#include <gsl/gsl_sf_bessel.h>
#include "../auxiliary/bessel_deriv_zero.hpp"
#include "../auxiliary/read_data.hpp"

#include <vector>
#include <assert.h>

//#include "../complex_matrix.hpp"


template<class REAL>
class FBDecompositionSimplified
{
public:
	typedef std::size_t size_type;
	typedef REAL number_type;

	//constructor
	FBDecompositionSimplified(
	  Model<number_type> model
	, number_type rMax
	 )
	: model_(model)
	, rMax_(rMax)
	, pi_(3.14159265359)
	, normalization_(1.)
	, Nm_(128)
	, Nr_(50)
	, b_percentiles_edges_(0)
	, b_percentiles_prob_(0)
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
		number_type x1;
		number_type y1;
		number_type angle = 0;
		x1 = x*cos(angle) + y*sin(angle);
		y1 = (-1.)*x*sin(angle) + y*cos(angle);
		x = x1;
		y = y1;

		return model_.OnePoint(x, y)/normalization_;
	}

	// Read in data on impact parameter distribution
	void get_impact_parameter_distribution(std::string filename)
	{
		std::vector<std::vector<number_type>> input(2);
		read_data(filename, input);

		b_percentiles_edges_ = input[0];
		b_percentiles_prob_ = input[1];
	}

	// Getter for two-Point correlation function
	number_type TwoPoint(number_type r, number_type b, number_type dPhi)
	{
		return model_.TwoPoint(r, b, dPhi);
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
		// Compute certain grid points of rho and generate a spline
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

		if (m == 0)
		{
			if (l == 1)
			{
				return 0;
			}
			else
			{
				return gsl_matrix_get(bessel_deriv_zeros_, m, l-2);	
			}
		}
		else
		{
			return gsl_matrix_get(bessel_deriv_zeros_, m, l-1);	
		}

	}

	// Evaluate appropriate Bessel function
	number_type Bessel(int m, size_type l, number_type r)
	{
		return gsl_sf_bessel_Jn(m, Bessel_zero(m, l)*rho(r));
	}



	// evaluate Bessel coefficient c_ml (if one uses J' zeros)
	number_type c (int m, size_type l)
	{
		number_type zero = Bessel_zero(m, l);
		number_type result;
		if ((m==0) && (l==1))
		{
			result = 1./2.;	
		}
		else
		{
			result = gsl_sf_bessel_Jn(m, zero)*gsl_sf_bessel_Jn(m, zero)*(zero*zero - m*m)/2./zero/zero;
		}
		return result; 
	}


	
	// my own integration routine for the r1, r2 integration
	number_type integ_r(int m1, int l1, int m2, int l2, number_type phi, size_type centrality_min)
	{

		// do integral over impact param by two evaluations

		number_type result_impact = 0;

		for (size_type i = 0; i < 2; ++i) // integral over impact param
		{
			number_type b = b_percentiles_edges_[centrality_min+i];
			number_type p = b_percentiles_prob_[centrality_min+i];

			number_type integral_r = 0;
			
			// compute integral over r at fixed phi
			// using trapezoidal rule
			size_type N = 25*2;
			number_type width = (rMax_-0.2)/N;
			for (size_type j = 0; j <= N; ++j) // loop over r
			{
				number_type r = (rMax_-0.2)*j/N;
				number_type f = p*r*Bessel(m1, l1, r)*Bessel(m2, l2, r)*TwoPoint(r, b, phi);
	
				if ((j == 0) || (j == N))
				{
					integral_r += f/2;
				}
				else
				{
					integral_r += f;
				}
			}
			integral_r *= width;

			result_impact += integral_r/2; // trapezoidal rule with N=2
		}

		result_impact *= (b_percentiles_edges_[centrality_min+1]-b_percentiles_edges_[centrality_min]);

		result_impact *= 100.; // since we are integrating the prob density over a percentile
	

		return result_impact;
	}

	// compute structure coefficient B_{l1l2}^{m1,m2, m } as a function of m
	void get_B(int l1, int l2, int m1, int m2, size_type centrality_min, std::vector<number_type> & Bm)
	{
		assert(Bm.size() == 0);
		get_Bessel_deriv_zeros(std::max(abs(m1), abs(m2))+1, std::max(l1, l2));
		//Compute FFT of a function E_l1l2^(m1m2)(phi1, phi2)
		
		// Step 1: Compute FFT with respect to phi
		
	
		// Do FFT with respect to phi
		number_type* fft;
		fft = new number_type[Nm_];

		for (size_type i = 0; i < Nm_; ++i)
		{
			number_type phi = 2.* pi_ * i / Nm_;
			fft[i] = integ_r(m1, l1, m2, l2, phi, centrality_min);
		}

		gsl_fft_real_radix2_transform(fft, 1, Nm_);
		
		// save result
		for (size_type i = 0; i < Nm_; ++i)
		{
			Bm.push_back(fft[i]/Nm_);
		}


		// Step 2: Scale result
		for (size_type i = 0; i < Bm.size(); ++i)
		{
			Bm[i] *= 1./2/pi_/c(m1, l1)/c(m2, l2);
		}

	}

	number_type get_IPSM(int l, int m)
	{
		size_type N = 170;
		number_type result = 1./2/pi_/pi_/N/c(m, l);
		if (m%2 == 1)
		{
			result = -result;
		}
		return result;
	}

	// compute two-mode correlator <e_l1^(m1) e_l2^(m2)>
	number_type TwoMode(int m1, int l1, int m2, int l2, size_type centrality_min)
	{
		get_Bessel_deriv_zeros(std::max(abs(m1), abs(m2))+1, std::max(l1, l2));
		
		number_type result;
			
		std::vector<number_type> Bm(0);
		get_B(l1, l2, m1, m2, centrality_min, Bm);

		// return real part of the result
		int m = m1+m2;
		if (m == 0)
		{
			result = Bm[0];
		}
		else if (m > 0)
		{
			result = Bm[m];
		}
		else
		{
			result = Bm[-m];
		}
		
		return result;

	}


	
	

	number_type pi_;
public:
	Model<number_type> model_;

private:
	number_type rMax_;
	size_type Nr_;
	size_type Nm_;
	number_type normalization_;

	const gsl_interp_type* W_interpolator_;
	gsl_spline* W_spline_;
	gsl_interp_accel* W_acc_;

	const gsl_interp_type* rho_interpolator_;
	gsl_spline* rho_spline_;
	gsl_interp_accel* rho_acc_;

	// vector of percentile edges for impact parameter
	std::vector<number_type> b_percentiles_edges_;
	// corresponding probability density
	std::vector<number_type> b_percentiles_prob_;

	gsl_matrix* bessel_deriv_zeros_;

};

















// struct to feed radius to phi integration function
template<typename number_type>
struct integ_one_phi_params
{
	number_type radius;
	FBDecompositionSimplified<number_type>* pt_FBDecomposition;
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
number_type FBDecompositionSimplified<number_type>::integ_one_phi(number_type r, number_type phi_lower, number_type phi_upper)
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
	FBDecompositionSimplified<number_type>* pt_FBDecomposition;
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
number_type FBDecompositionSimplified<number_type>::integ_one_r(number_type r_lower, number_type r_upper)
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
	FBDecompositionSimplified<number_type>* pt_FBDecomposition;
};

// gsl function to call to integrate over W(r)r
double integ_W_r_wrapper(double r, void* params)
{
	integ_W_r_params<double>* params_ = (integ_W_r_params<double>*) params;

	return (params_->pt_FBDecomposition->W(r))*r;
}

// integrate W(r)r over r
template<class number_type>
number_type FBDecompositionSimplified<number_type>::integ_W_r(number_type r_lower, number_type r_upper)
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

// // struct for integ_r2
// template<typename size_type, typename number_type>
// struct integ_r2_params
// {
// 	size_type m2;
// 	size_type l2;
// 	number_type phi1;
// 	number_type phi2;
// 	number_type r1;
// 	FBDecomposition<number_type>* pt_FBDecomposition;
// };

// // gsl function to call to integrate over r2
// double integ_r2_wrapper(double r, void* params)
// {
// 	typedef int size_type;
// 	typedef double number_type;
// 	integ_r2_params<size_type, number_type>* params_ = (integ_r2_params<size_type, number_type>*) params;
// 	size_type m2_ = params_->m2;
// 	size_type l2_ = params_->l2;
// 	number_type phi1_ = params_->phi1;
// 	number_type phi2_ = params_->phi2;
// 	number_type r1_ = params_->r1;
// 	return (params_->pt_FBDecomposition->weight(m2_, l2_, r))*(params_->pt_FBDecomposition->TwoPoint(r1_, phi1_, r, phi2_));
// }


// // integrate two-point function over r2 at fixed phi1, phi2, r1, m2, l2
// // and Besselfunction in r2
// template<class number_type>
// number_type FBDecomposition<number_type>::integ_r2(int m2, int l2, number_type phi1, number_type phi2, number_type r1)
// {
// 	typedef int size_type;
// 	gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);
// 	number_type integral;
// 	number_type error;
// 	integ_r2_params<size_type, number_type> params = {m2, l2, phi1, phi2, r1, this};

// 	gsl_function F;
// 	F.function = &integ_r2_wrapper;
// 	F.params = &params;
// 	gsl_integration_qags(&F, 0, (this->get_rMax()-0.2), 1e-5, 1e-5, 1000, w, &integral, &error);
// 	//std::cout << " Radius: " << r << " Index: " << k <<  " Integral: " << integral << " Error abs: " << error << " Error rel: " << error/integral << "\n"; 
// 	gsl_integration_workspace_free(w);
// 	return integral;
// }

// // struct for integ_r1
// template<typename size_type, typename number_type>
// struct integ_r1_params
// {
// 	size_type m1;
// 	size_type l1;
// 	size_type m2;
// 	size_type l2;
// 	number_type phi1;
// 	number_type phi2;
// 	FBDecomposition<number_type>* pt_FBDecomposition;
// };

// // gsl function to call to integrate over r1
// double integ_r1_wrapper(double r, void* params)
// {
// 	typedef int size_type;
// 	typedef double number_type;
// 	integ_r1_params<size_type, number_type>* params_ = (integ_r1_params<size_type, number_type>*) params;
// 	size_type m1_ = params_->m1;
// 	size_type l1_ = params_->l1;
// 	size_type m2_ = params_->m2;
// 	size_type l2_ = params_->l2;
// 	number_type phi1_ = params_->phi1;
// 	number_type phi2_ = params_->phi2;
// 	return (params_->pt_FBDecomposition->weight(m1_, l1_, r))*(params_->pt_FBDecomposition->integ_r2(m2_, l2_, phi1_, phi2_, r));
// }


// // integrate two-point function over r1 at fixed phi1, phi2, m2, l2
// // and Besselfunction in r1
// template<class number_type>
// number_type FBDecomposition<number_type>::integ_r1(int m1, int l1, int m2, int l2, number_type phi1, number_type phi2)
// {
// 	typedef int size_type;
// 	gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);
// 	number_type integral;
// 	number_type error;
// 	integ_r1_params<size_type, number_type> params = {m1, l1, m2, l2, phi1, phi2, this};

// 	gsl_function F;
// 	F.function = &integ_r1_wrapper;
// 	F.params = &params;
// 	gsl_integration_qags(&F, 0, (this->get_rMax()-0.2), 1e-5, 1e-5, 1000, w, &integral, &error);
// 	//std::cout << " Radius: " << r << " Index: " << k <<  " Integral: " << integral << " Error abs: " << error << " Error rel: " << error/integral << "\n"; 
// 	gsl_integration_workspace_free(w);
// 	return integral;
// }

#endif