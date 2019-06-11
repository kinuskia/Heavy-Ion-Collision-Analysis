#ifndef COLLISION_HPP
#define COLLISION_HPP

#include <vector>
#include "read_data.hpp"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_fft_real.h>
#include <algorithm>
#include "to_file.hpp"
#include "statistics.hpp"
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sf_bessel.h>
#include "complex_matrix.hpp"
#include "bessel_deriv_zero.hpp"
#include <fstream>
#include <random>
#include <cmath>


// Functions to transform polar coordinates to cartesian coordinates
template<typename number_type>
number_type to_x(number_type r, number_type phi)
{
	return r*cos(phi);
}

template<typename number_type>
number_type to_y(number_type r, number_type phi)
{
	return r*sin(phi);
}


template <class REAL>
class Collision
{
public:
	typedef std::size_t size_type;
	typedef REAL number_type;
	//constructor
	Collision(number_type grid_max, number_type grid_step)
	: grid_max_(grid_max)
	, grid_step_(grid_step)
	, grid_size_(size_type(2.*grid_max/grid_step))
	, profiles_(0)
	, random_reaction_plane_(false)
	, angles_(0)
	, m_profiles_(0)
	, z_profiles_(0)
	, impact_parameters_(0)
	, n_participants_(0)
	, multiplicities_(0)
	, percentiles_(0)
	, percent_mult_(0)
	, interpolators_(0)
	, splines_(0)
	, xacc_(0)
	, yacc_(0)
	, r_interpolators_(0)
	, r_splines_(0)
	, racc_(0) 
	, Nm_(128)
	, Nr_(size_type(grid_max/grid_step))
	, W_interpolators_(0)
	, W_splines_(0)
	, W_acc_(0)
	, normalizations_(0)
	{
		// initialize grid sites
		x_sites_ = new number_type[grid_size_];
		y_sites_ = new number_type[grid_size_];
		r_sites_ = new number_type[Nr_];
		R_ = grid_max_ - grid_step_;
		for (size_type i = 0; i < grid_size_; ++i)
		{
			x_sites_[i] = x_from_index(i);
			y_sites_[i] = -y_from_index(i);
			r_sites_[i] = R_ * i / Nr_;
		}
	}

	// getter functions
	size_type Nm() const
	{
		return Nm_;
	}

	// free profiles
	void free()
	{
		for (size_type i = 0; i < profiles_.size(); ++i)
		{
			gsl_matrix_free(profiles_[i]);
			gsl_matrix_free(m_profiles_[i]);
		}
	}

	// function to generate the filename of data file k
	std::string get_filename(std::string filename_begin, std::string fileformat, size_type k, size_type n_files)
	{
		std::string filename = filename_begin;
		
		// Step 1: Determine the number of digits to display string number
		size_type n_digits = log(1.0*(n_files-1))/log(10.) + 1;
		
		// Step 2: Determine the number of leading zeros needed
		size_type n_digits_index;
		if (k == 0)
		{
			n_digits_index = 1;
		}
		else
		{
			n_digits_index = log(1.0*k)/log(10.) + 1;
		}
		size_type leading_zeros = n_digits - n_digits_index;
		
		// Step 3: Create file name
		for (size_type i = 0; i < leading_zeros; ++i)
		{
			filename += "0";
		}
		filename += std::to_string(k);
		filename += fileformat;
		return filename;
	}


	// // read in single file
	// void read_in(std::string filename)
	// {
	// 	gsl_matrix* data = gsl_matrix_alloc(grid_size_, grid_size_);
	// 	read_data(filename, data, 8);
	// 	profiles_.push_back(data);
	// }

	// read in data files
	void read_in(std::string filename_begin, std::string fileformat, size_type n_files, number_type exponent, bool import_profiles = true)
	{
		// initialize random devices to sample random reaction plane angles
		std::random_device rd;
		std::mt19937 gen(rd());
		std::uniform_real_distribution<> dis(0, 2*3.1415926);

		for (size_type k = 0; k < n_files; ++k)
		{
			std::string filename = get_filename(filename_begin, fileformat, k, n_files);
			number_type impact_param;
			number_type n_part;
			number_type mult;
			if (import_profiles)
			{	
				gsl_matrix* data = gsl_matrix_alloc(grid_size_, grid_size_);
				read_data(filename, data, impact_param, n_part, mult);
				// compute the power of the profiles
				for (size_type i = 0; i < data->size1; ++i)
				{
					for (size_type j = 0; j < data->size2; ++j)
					{
						number_type density = gsl_matrix_get(data, i, j);
						density = pow(density, exponent);
						gsl_matrix_set(data, i, j, density);	
					}
				}
				profiles_.push_back(data);
			}
			else
			{
				read_data(filename, impact_param, n_part, mult);	
			}
			
			impact_parameters_.push_back(impact_param);
			n_participants_.push_back(n_part);
			multiplicities_.push_back(mult);

			if (random_reaction_plane_)
			{
				number_type angle;
				angle = dis(gen);
				angles_.push_back(angle);
			}
		}
	}

	// Save header data of trento files in a separate text file

	void collision_specs_to_file(std::string filename) const
	{
		std::vector<std::vector<number_type>> columns(3);
		columns[0] = impact_parameters_;
		columns[1] = n_participants_;
		columns[2] = multiplicities_;
		to_file(filename, columns);
	}

	// Overload for specific centrality class
	// Save header data of trento files in a separate text file

	void collision_specs_to_file(std::string filename, size_type centrality_index ) const
	{
		assert(centrality_index > 0 && centrality_index < percentiles_.size());
		std::vector<std::vector<number_type>> columns(3);
		// columns[0] : impact_parameters
		// columns[1] : n_participants
		// columns[2] : multiplicities
		for (size_type k = 0; k < impact_parameters_.size(); ++k)
		{
			if (is_in_centrality_class(k, centrality_index))
			{
				columns[0].push_back(impact_parameters_[k]);
				columns[1].push_back(n_participants_[k]);
				columns[2].push_back(multiplicities_[k]);
			}
			else
			{
				continue;
			}
		}


		to_file(filename, columns);
	}

	// compute histogram data points of multiplicity samples
	void histogram_mult(std::string filename, size_type N_bins, bool normed)
	{
		// create edges and bin vector
		std::vector<std::size_t> frequencies(N_bins, 0);
		std::vector<number_type> samples = multiplicities_;
		number_type min = *std::min_element(samples.begin(), samples.end());
		number_type max = *std::max_element(samples.begin(), samples.end());

		std::vector<number_type> edges(N_bins+1);
		for (size_type i = 0; i < edges.size(); ++i)
		{
			edges[i] = min + (max-min)/N_bins*i;
		}

		// compute histogram
		histogram(samples, edges, frequencies);

		number_type integral = 0;
		for (size_type i = 0; i < N_bins; ++i)
		{
			number_type width = edges[i+1]-edges[i];
			integral += frequencies[i]*width;
		}

		if (!normed)
		{
			integral = 1.0;
		}

		// save to text file
		std::vector<number_type> x(N_bins);
		std::vector<number_type> y(N_bins);
		std::vector<number_type> dy(N_bins);
		for (size_type i = 0; i < x.size(); ++i)
		{
			x[i] = (edges[i+1]+edges[i])/2;
			y[i] = number_type(frequencies[i])/integral;
			dy[i] = sqrt(number_type(frequencies[i]))/integral;
		}

		std::vector<std::vector<number_type>> columns(3);
		columns[0] = x;
		columns[1] = y;
		columns[2] = dy;

		to_file(filename, columns);
	}

	// compute histogram data points of impact parameters
	void histogram_b(std::string filename, size_type N_bins, bool normed)
	{
		// create edges and bin vector
		std::vector<std::size_t> frequencies(N_bins, 0);
		std::vector<number_type> samples = impact_parameters_;
		number_type min = *std::min_element(samples.begin(), samples.end());
		number_type max = *std::max_element(samples.begin(), samples.end());

		std::vector<number_type> edges(N_bins+1);
		for (size_type i = 0; i < edges.size(); ++i)
		{
			edges[i] = min + (max-min)/N_bins*i;
		}

		// compute histogram
		histogram(samples, edges, frequencies);

		number_type integral = 0;
		for (size_type i = 0; i < N_bins; ++i)
		{
			number_type width = edges[i+1]-edges[i];
			integral += frequencies[i]*width;
		}

		if (!normed)
		{
			integral = 1.0;
		}

		// save to text file
		std::vector<number_type> x(N_bins);
		std::vector<number_type> y(N_bins);
		std::vector<number_type> dy(N_bins);
		for (size_type i = 0; i < x.size(); ++i)
		{
			x[i] = (edges[i+1]+edges[i])/2;
			y[i] = number_type(frequencies[i])/integral;
			dy[i] = sqrt(number_type(frequencies[i]))/integral;
		}

		std::vector<std::vector<number_type>> columns(3);
		columns[0] = x;
		columns[1] = y;
		columns[2] = dy;

		to_file(filename, columns);
	}


	// compute certain multiplicity percentiles
	void get_percentiles(const std::vector<number_type> & percentiles)
	{
		std::vector<number_type> multip_copy = multiplicities_;
		size_type N = multiplicities_.size();
		for (size_type i = 0; i< percentiles.size(); ++i)
		{
			normalizations_.push_back(1.);
			percentiles_.push_back(percentiles[i]);
			size_type dn = (N-1)*(100-percentiles[i])/100;
			std::nth_element(multip_copy.begin(), multip_copy.begin()+dn, multip_copy.end());
			percent_mult_.push_back(multip_copy[dn]);
			std::cout << percentiles_[i] << " " << percent_mult_[i] << "\n";
		}
	}

	void random_reaction_plane(bool yes)
	{
		if (yes)
		{
			random_reaction_plane_ = true;
		}
		else
		{
			random_reaction_plane_ = false;
		}
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

	// print Bessel derivative zeros to text file
	void print_Bessel_deriv_zeros(size_type m, size_type l, std::string filename)
	{
		get_Bessel_deriv_zeros(m, l);
		to_file(filename, bessel_deriv_zeros_);
	}

	number_type Bessel_zero(size_type m, size_type l)
	{
		return gsl_matrix_get(bessel_deriv_zeros_, m, l-1);
		//return gsl_sf_bessel_zero_Jnu(m, l);
	}


	// initialize xy-interpolation
	void initialize_xy_interpolation(const gsl_interp2d_type* interp_method)
	{
		for (size_type k = 0; k < profiles_.size(); ++k)
		{
			const gsl_interp2d_type* interpolator = interp_method;
			interpolators_.push_back(interpolator);
			number_type* za = new number_type[grid_size_*grid_size_];
			z_profiles_.push_back(za);
			gsl_spline2d* spline = gsl_spline2d_alloc(interpolators_[k], grid_size_, grid_size_);
			splines_.push_back(spline);
			gsl_interp_accel* xacc = gsl_interp_accel_alloc();
			gsl_interp_accel* yacc = gsl_interp_accel_alloc();
			xacc_.push_back(xacc);
			yacc_.push_back(yacc);
			// set interpolation grid values
			for (size_type i = 0; i < grid_size_; ++i)
			{
				for (size_type j = 0; j < grid_size_; ++j)
				{
					gsl_spline2d_set(splines_[k], z_profiles_[k], i, j, gsl_matrix_get(profiles_[k], grid_size_ -1 - i, j));
				}
			}
			gsl_spline2d_init(splines_[k], x_sites_, y_sites_, z_profiles_[k], grid_size_, grid_size_);

			// free profile
			gsl_matrix_free(profiles_[k]);
		}

	}

	// initialize r-interpolation
	/*
		The interpolation objects of the m'th mode of the k'th profile are saved at the index:
		index = k*Nm_ + m
		-> conversion:
		k = index/Nm_
		m = index%Nm_
	*/
	void initialize_r_interpolation(const gsl_interp_type* interp_method)
	{
		size_type N = profiles_.size();
		for (size_type k = 0; k < N; ++k)
		{
			for (size_type m = 0; m < Nm_; ++m)
			{
				size_type index = k * Nm_ + m; // last index value: (N-1)*Nm_ + Nm_-1 = N*Nm_-1
				const gsl_interp_type* interpolator = interp_method;
				r_interpolators_.push_back(interpolator);
				gsl_spline* spline = gsl_spline_alloc(r_interpolators_[index], Nr_);
				r_splines_.push_back(spline);
				gsl_interp_accel* racc = gsl_interp_accel_alloc();
				racc_.push_back(racc);

				// set interpolation values
				number_type* e_sites;
				e_sites = new number_type[Nr_];
				for (size_type i = 0; i < Nr_; ++i)
				{
					e_sites[i] = gsl_matrix_get(m_profiles_[k], m, i);
				} 


				// initialize spline object
				gsl_spline_init(r_splines_[index], r_sites_, e_sites, Nr_);
			}

			// free m_profiles
			gsl_matrix_free(m_profiles_[k]);
		}

	}

	// normalize matrix data (integral = 1)
	void normalize(number_type energy)
	{
		for (size_type k = 0; k < profiles_.size(); ++k)
		{
			// compute total energy
			number_type sum = 0;
			for (size_type i = 0; i < profiles_[k]->size1; ++i)
			{
				for (size_type j = 0; j < profiles_[k]->size2; ++j)
				{
					sum += gsl_matrix_get(profiles_[k], i, j);
				}
			}
			sum *= grid_step_*grid_step_;

			std::cout << sum << "\n";
			
			// rescale
			gsl_matrix_scale(profiles_[k], energy/sum);
		}
	}

	// Integrate profile k
	number_type integrated_density(size_type k)
	{
		// compute total energy
		number_type sum = 0;
		for (size_type i = 0; i < profiles_[k]->size1; ++i)
		{
			for (size_type j = 0; j < profiles_[k]->size2; ++j)
			{
				sum += gsl_matrix_get(profiles_[k], i, j);
			}
		}
		sum *= grid_step_*grid_step_;

		return sum;
		
	}

	// centralize matrix data (barycenter at N/2,N/2)
	void centralize()
	{
		for (size_type k = 0; k < profiles_.size(); ++k)
		{
			// compute barycentre
			number_type x_bary = 0;
			number_type y_bary = 0;
			number_type integral = integrated_density(k);
			for (size_type i = 0; i < profiles_[k]->size1; ++i)
			{
				for (size_type j = 0; j < profiles_[k]->size2; ++j)
				{
					number_type x = j;
					number_type y = i;
					x_bary += gsl_matrix_get(profiles_[k], i, j)*x/integral;
					y_bary += gsl_matrix_get(profiles_[k], i, j)*y/integral;
				}
			}
			x_bary *= grid_step_*grid_step_;
			y_bary *= grid_step_*grid_step_;


			// make temporary copy
			gsl_matrix* m_copy = gsl_matrix_alloc(profiles_[k]->size1, profiles_[k]->size2);
			gsl_matrix_memcpy(m_copy, profiles_[k]);

			// center the barycentre
			int dx = round((1.0*m_copy->size2-1)/2-x_bary);
			int dy = round((1.0*m_copy->size1-1)/2-y_bary);
			//std::cout << "x: " << xB << " y: " << yB << "\n";
			
			for (int i = 0; i < m_copy->size1; ++i)
			{
				for (int j = 0; j < m_copy->size2; ++j)
				{
					// deal with out-of bound indices
					if (i-dy < 0 || j-dx < 0 || i-dy >= m_copy->size2 || j-dx >= m_copy->size1)
					{
						gsl_matrix_set(profiles_[k], i, j, 0.);
					}
					else
					{
						gsl_matrix_set(profiles_[k], i, j, gsl_matrix_get(m_copy, i-dy, j-dx));
					}
				}
			}
			gsl_matrix_free(m_copy);
		}
	}

	// Get normalization constants for each centrality class
	void getNormalizations()
	{
		for (size_type c = 1; c < percentiles_.size(); ++c)
		{
			std::cout << normalizations_[c-1] << "\n";
			normalizations_[c-1] = (2.*integ_r_Wr(c, 0, R_-grid_step_));
			std::cout << normalizations_[c-1] << "\n";
		}

		// normalize m_profiles
		for (size_type c = 1; c < percentiles_.size(); ++c)
		{
			for (size_type k = 0; k < m_profiles_.size(); ++k)
			{
				// Skip profile if not in centrality class
				if (!is_in_centrality_class(k, c))
				{
					continue;
				}

				gsl_matrix_scale(m_profiles_[k], 1./normalizations_[c-1]);

			}
		}
	}


	// get reaction plane angle
	void getReactionPlane(std::string filename)
	{
		std::vector<number_type> angles(0);

		for (size_type k = 0; k < profiles_.size(); ++k)
		{
			// compute entries of inertia tensor
			number_type Qxx = 0;
			number_type Qxy = 0;
			number_type Qyy = 0;


			for (size_type i = 0; i < profiles_[k]->size1; ++i)
			{
				for (size_type j = 0; j < profiles_[k]->size2; ++j)
				{
					number_type x = x_from_index(j);
					number_type y = y_from_index(i);
					Qxx += gsl_matrix_get(profiles_[k], i, j)*(y*y);
					Qxy += gsl_matrix_get(profiles_[k], i, j)*(-x*y);
					Qyy += gsl_matrix_get(profiles_[k], i, j)*(x*x);
				}
			}
			Qxx *= grid_step_*grid_step_*grid_step_;
			Qxy *= grid_step_*grid_step_*grid_step_;
			Qyy *= grid_step_*grid_step_*grid_step_;


			// compute angle by which the profile needs to be rotated clockwise
			number_type val = (Qyy-Qxx)/2./Qxy;
			number_type angle = atan(val + sqrt(1+val*val)*sign(Qxy));

			angles.push_back(angle);
		}
		std::vector<std::vector<number_type>> angle_column(1);
		angle_column[0] = angles;
		to_file(filename, angle_column);
	}

	// print ensemble averaged profiles
	void print_averaged_profiles(std::string filename)
	{
		for (size_type c = 1; c < percentiles_.size(); ++c) // loop over centrality classes
		{
			gsl_matrix* sum = gsl_matrix_alloc(profiles_[0]->size1, profiles_[0]->size2);
			gsl_matrix_scale(sum, 0.0);
			size_type counter = 0;
			for (size_type k = 0; k < profiles_.size(); ++k)
			{
				if (is_in_centrality_class(k, c))
				{
					// add profile data to sum
					gsl_matrix_add(sum, profiles_[k]);
					counter++;
				}
				else
				{
					continue;
				}
			}
			// divide by the number of elements to obtain average
			gsl_matrix_scale(sum, number_type(1./counter));

			// save ensemble average in text file
			std::string outfilename = filename;
			outfilename += "_" + std::to_string(size_type(percentiles_[c-1])) + "-" + std::to_string(size_type(percentiles_[c]));
			outfilename += ".txt";
			to_file(outfilename, sum);
			gsl_matrix_free(sum);
		}

	}

	// Functions to convert index pair (i,j) to cartesian coordinates (x,y)
	number_type x_from_index(size_type j) const
	{
		return -grid_max_+grid_step_/2 + grid_step_*j;
	}

	number_type y_from_index(size_type i) const
	{
		return grid_max_-grid_step_/2 - grid_step_*i;
	}


	// integrate over phi at fixed r with gsl integration routine
	// (only declaration, it is defined outside this class, see below)
	number_type integ_phi(size_type k, number_type r, number_type phi_lower, number_type phi_upper);

	// integrate over r at fixed Fourier mode m and Bessel(m,l) as weight function with gsl integration routine
	// (only declaration, it is defined outside this class, see below)
	number_type integ_fourier_r(size_type k, size_type centrality_index, size_type m, size_type l, number_type r_lower, number_type r_upper);
	

	// integrate W(r) r epsilon(r, phi) over r for given profile
	number_type integ_r(size_type k, number_type r_lower, number_type r_upper);

	// integrate W(r) r over r for given centrality index
	number_type integ_r_Wr(size_type centrality_index, number_type r_lower, number_type r_upper);

	// interpolation of profil k
	number_type interpolate(size_type k, number_type x, number_type y) 
	{
		// consider random reaction plane angle
		if (random_reaction_plane_)
		{
			number_type x1;
			number_type y1;
			number_type angle = angles_[k];
			x1 = x*cos(angle) + y*sin(angle);
			y1 = (-1.)*x*sin(angle) + y*cos(angle);
			x = x1;
			y = y1;
		}
		return gsl_spline2d_eval(splines_[k], x, y, xacc_[k], yacc_[k]);
	}

	// r-interpolation of m-th Fourier mode of profile k
	number_type interpolate_fourier(size_type k, size_type m, number_type r)
	{
		size_type index = k * Nm_ + m;
		return gsl_spline_eval(r_splines_[index], r, racc_[index]);
	}

	// interpolate weight function
	number_type W(number_type r, size_type centrality_index)
	{
		return gsl_spline_eval(W_splines_[centrality_index-1], r, W_acc_[centrality_index-1])/normalizations_[centrality_index-1];
	}

	// Evaluate appropriate Bessel function
	number_type Bessel(size_type m, size_type l, number_type r, size_type centrality_index)
	{
		// // Verify the zero crossing needed is available
		// assert(l>0);
		// assert(l <= bessel_deriv_zeros_->size2);
		// assert(m < bessel_deriv_zeros_->size1);

		return gsl_sf_bessel_Jn(m, Bessel_zero(m, l)*rho(r, centrality_index));
	}

	// evaluate Bessel coefficient c_ml (if one uses J' zeros)
	number_type c (size_type m, size_type l)
	{
		number_type zero = Bessel_zero(m, l);
		return gsl_sf_bessel_Jn(m, zero)*gsl_sf_bessel_Jn(m, zero)*(zero*zero - m*m)/2./zero/zero; 

		// if (m == 0 && l == 1)
		// {
		// 	return 0.5;
		// }
		// else if (m == 0 && l > 1)
		// {
		// 	number_type zero = Bessel_zero(m, l-1);
		// 	return gsl_sf_bessel_Jn(m, zero)*gsl_sf_bessel_Jn(m, zero)*(zero*zero - m*m)/2./zero/zero;
		// }
		// else
		// {
		// 	number_type zero = Bessel_zero(m, l);
		// 	return gsl_sf_bessel_Jn(m, zero)*gsl_sf_bessel_Jn(m, zero)*(zero*zero - m*m)/2./zero/zero;
		// }
		
	}

	// // evaluate Bessel coefficient c_ml (if one uses J zeros)
	// number_type c (size_type m, size_type l)
	// {
	// 	number_type zero = Bessel_zero(m, l); 
	// 	// value of J' at zero
	// 	number_type h = 0.001;
	// 	number_type deriv = (gsl_sf_bessel_Jn(m, zero+h)-gsl_sf_bessel_Jn(m, zero-h))/2./h;
	// 	return (deriv*deriv)/2.;
	// }


	// map r to normalized variable rho \in (0, 1) using W(r)
	// only declared, definition outside class
	number_type rho(number_type r, size_type centrality_index);

	
	// declaration for test function, see below
	number_type integ_test();


	// average over azimuthal angle save to data storage
	void average_azimuthal(gsl_vector* radii, gsl_vector* mean_vec, gsl_vector* mean_error)
	{
		size_type N = radii->size;

		// For each profile compute e_average as function of r
		for (size_type l = 0; l < N; ++l) // loop over radii
		{
			number_type r = number_type(l)/N*R_;
			gsl_vector* e_ensemble = gsl_vector_alloc(profiles_.size());
			for (size_type k = 0; k < profiles_.size(); ++k) // loop over profiles
			{
				if (k == 0)
				{
					gsl_vector_set(radii, l, r);
				}
				
				number_type e_average = integ_phi(k, r, 0., 2.*3.1415926)/(2.*3.1415926);

				gsl_vector_set(e_ensemble, k, e_average);
			}
			// Compute mean and uncertainty of the fluctuations
			number_type e_mean;
			number_type e_err;
			mean(e_ensemble, e_mean, e_err);
			gsl_vector_set(mean_vec, l, e_mean);
			gsl_vector_set(mean_error, l, e_err);

			gsl_vector_free(e_ensemble);
		}
	}

	// // Compute ensemble-averaged expectation value of background coefficient
	// void getOnePointFunction_00(size_type centrality_index, number_type & coeff_mean, number_type & coeff_err, std::time_t start)
	// {

	// 	std::vector<number_type> coeffs(0);

	// 	for (size_type k = 0; k < profiles_.size(); ++k)
	// 	{
	// 		// Skip profile if not in centrality class
	// 		if (!is_in_centrality_class(k, centrality_index))
	// 		{
	// 			continue;
	// 		}
	// 		number_type current_coeff;
	// 		current_coeff = getBackgroundCoefficient(k);
	// 		coeffs.push_back(current_coeff);
	// 	}

	// 	mean(coeffs, coeff_mean, coeff_err);
	// 	std::cout << "class size: " << coeffs.size() << "\n";
	// 	std::time_t current_time = std::time(nullptr);
	// 	std::cout << current_time-start << "s: " << "Backgound one-mode-fluctuation has been computed.\n";
	// }

	// Compute one-point correlation function (=expectation value e_ml) for given centrality class
	void getOnePointFunction(size_type centrality_index, complex_matrix<number_type> & coeffs_mean, complex_matrix<number_type> & coeffs_err, std::time_t start)
	{


		// // Bessel decomposition

		size_type mMax = coeffs_mean.rowsize()-1;
		size_type lMax = coeffs_mean.colsize();

		// computed necessary bessel zeros
		get_Bessel_deriv_zeros(mMax+1, lMax);

		// compute Fourier-Bessel coeffs

		std::vector<complex_matrix<number_type>> Fourier_Bessel_coeffs(0);

		for (size_type k = 0; k < profiles_.size(); ++k)
		{
			// Skip profile if not in centrality class
			if (!is_in_centrality_class(k, centrality_index))
			{
				continue;
			}
			complex_matrix<number_type> current_coeff(mMax+1, lMax);
			BesselDecomposeProfile(k, centrality_index, current_coeff);

			
			Fourier_Bessel_coeffs.push_back(current_coeff);
		}

		std::time_t current_time = std::time(nullptr);
		std::cout << current_time-start << "s: " << "Profiles have been Bessel-decomposed.\n"; 

		// Compute ensemble mean of Fourier-Bessel coefficients
		mean(Fourier_Bessel_coeffs, coeffs_mean, coeffs_err);

		// CAREFUL: commented lines no longer work (is now centrality class specific)!!

		// // print energy-r curves for some numbers m
		// for (size_type m = 0; m < 4; ++m)
		// {
		// 	gsl_vector* radii = gsl_vector_alloc(Nr_);
		// 	gsl_vector* e_mean = gsl_vector_alloc(Nr_);
		// 	gsl_vector* e_err = gsl_vector_alloc(Nr_);

		// 	average_fourier(m, radii, e_mean, e_err);

		// 	std::string filename = "output/e_ave_fourier";
		// 	filename += std::to_string(m);
		// 	filename += ".txt";

		// 	to_file(filename, radii, e_mean, e_err);

		// 	gsl_vector_free(radii);
		// 	gsl_vector_free(e_mean);
		// 	gsl_vector_free(e_err);

		// }

		
		// // Print some bessel curves
		// size_type N = 100;
		// std::vector<number_type> radii(N);
		// std::vector<number_type> bessels11(N);
		// std::vector<number_type> bessels12(N);
		// std::vector<number_type> bessels13(N);
		// for (size_type i = 0; i < N; ++i)
		// {
		// 	radii[i] = (R_-grid_step_)*i/N;
		// 	bessels11[i] = Bessel(1, 1, radii[i])*W(radii[i])*R_*R_;
		// 	bessels12[i] = Bessel(1, 2, radii[i])*W(radii[i])*R_*R_;
		// 	bessels13[i] = Bessel(1, 3, radii[i])*W(radii[i])*R_*R_;
			
		// }
		// std::vector<std::vector<number_type>> data(4);
		// data[0] = radii;
		// data[1] = bessels11;
		// data[2] = bessels12;
		// data[3] = bessels13;
		// to_file("output/bessel1.txt", data);


	}

	// outsource certain evaluations of getOnePointFunction and getTwoPointFunction
	void initialize_n_point_evaluations(const gsl_interp_type* interpolation_method, std::time_t start)
	{
		/* Fourier-decompose the profiles */
		decompose_azimuthal();
		getNormalizations();
		initialize_r_interpolation(interpolation_method);


		std::time_t current_time = std::time(nullptr);
		std::cout << current_time-start << "s: " << "Profiles have been Fourier-decomposed.\n"; 	
	}

	// Method to compute a two-point correlation function of the form <e_ml eml'>
	void getTwoPointFunction_fixed_m(size_type m, size_type centrality_index, complex_matrix<number_type> & TwoPointFunction, complex_matrix<number_type> & TwoPointFunction_err, std::time_t start)
	{
		std::time_t current_time = std::time(nullptr);

		/* Bessel-decompose profile parts with quantum number m */

		size_type lMax = TwoPointFunction.colsize();

		// compute necessary bessel zeros
		get_Bessel_deriv_zeros(m+1, lMax);

		// compute two-point object e_ml e_ml' for each profile in the centrality class

		std::vector<complex_matrix<number_type>> Two_point_objects(0);

		for (size_type k = 0; k < profiles_.size(); ++k)
		{
			// Skip profile if not in centrality class
			if (!is_in_centrality_class(k, centrality_index))
			{
				continue;
			}


			complex_matrix<number_type> current_coeff(1, lMax);
			BesselDecomposeProfile_fixed_m(k, centrality_index, m, current_coeff);
			// fill current two-point object e_ml1 e*_ml2
			complex_matrix<number_type> current_two_point(lMax, lMax);
			for (size_type l1 = 0; l1 < current_two_point.rowsize(); ++l1)
			{
				for (size_type l2 = 0; l2 < current_two_point.colsize(); ++l2)
				{
					number_type eml1_real = current_coeff.get_real(0, l1);
					number_type eml1_imag = current_coeff.get_imag(0, l1);
					number_type eml2_real = current_coeff.get_real(0, l2);
					number_type eml2_imag = current_coeff.get_imag(0, l2)*(-1.); // minus sign because we want -m (complex conjugate) e_(-m,l) = (-1)^m e*_(m,l)

					if (m%2 == 1) // if m odd, we need a (-1)^m prefactor
					{
						eml2_real *= -1.;
						eml2_imag *= -1.;
					}
					
					number_type result_real = eml1_real*eml2_real - eml1_imag*eml2_imag;
					number_type result_imag = eml1_real*eml2_imag + eml2_real*eml1_imag;
					current_two_point.set_entry(l1, l2, result_real, result_imag);
				}
			}
			
			Two_point_objects.push_back(current_two_point);
		}

		current_time = std::time(nullptr);
		std::cout << current_time-start << "s: " << "Two-point objects have been computed.\n"; 

		// Compute ensemble mean of Fourier-Bessel coefficients
		mean(Two_point_objects, TwoPointFunction, TwoPointFunction_err);


	}

	// Method to compute a two-point correlation function of the form <e_(m1)l e(m2)l>
	void getTwoPointFunction_fixed_l(size_type l, size_type centrality_index, complex_matrix<number_type> & TwoPointFunction, complex_matrix<number_type> & TwoPointFunction_err, std::time_t start)
	{
		std::time_t current_time = std::time(nullptr);

		/* Bessel-decompose profile parts with quantum number m */

		size_type m_size = TwoPointFunction.rowsize();

		// compute necessary bessel zeros
		get_Bessel_deriv_zeros(m_size, l);

		// compute two-point object e_m1l e_m2l for each profile in the centrality class

		std::vector<complex_matrix<number_type>> Two_point_objects(0);

		for (size_type k = 0; k < profiles_.size(); ++k)
		{
			// Skip profile if not in centrality class
			if (!is_in_centrality_class(k, centrality_index))
			{
				continue;
			}


			complex_matrix<number_type> current_coeff(m_size, 1); 
			BesselDecomposeProfile_fixed_l(k, centrality_index, l, current_coeff);
			// fill current two-point object e_(m1)l e*_(-m2)l
			complex_matrix<number_type> current_two_point(m_size, 2*m_size-1); // will contain entries for m2 in (-mMax,...mMax)
			for (size_type m1 = 0; m1 < m_size; ++m1)
			{
				for (size_type m2 = 0; m2 < m_size; ++m2)
				{
					number_type em1l_real = current_coeff.get_real(m1, 0);
					number_type em1l_imag = current_coeff.get_imag(m1, 0);
					number_type em2l_real_pos = current_coeff.get_real(m2, 0);
					number_type em2l_imag_pos = current_coeff.get_imag(m2, 0); 
					number_type em2l_real_neg = em2l_real_pos;
					number_type em2l_imag_neg = em2l_imag_pos;

					// em2l looks different for negative m2
					em2l_imag_neg *= -1.; // minus sign because we want -m (complex conjugate) e_(-m,l) = (-1)^m e*_(m,l)
						
					if (m2%2 == 1) // if m odd, we need a (-1)^m prefactor
					{
						em2l_real_neg *= -1.;
						em2l_imag_neg *= -1.;
					}
					

					
					
					number_type result_real_pos = em1l_real*em2l_real_pos - em1l_imag*em2l_imag_pos;
					number_type result_imag_pos = em1l_real*em2l_imag_pos + em2l_real_pos*em1l_imag;
					number_type result_real_neg = em1l_real*em2l_real_neg - em1l_imag*em2l_imag_neg;
					number_type result_imag_neg = em1l_real*em2l_imag_neg + em2l_real_neg*em1l_imag;
					current_two_point.set_entry(m1,m_size-1+m2, result_real_pos, result_imag_pos);
					current_two_point.set_entry(m1,m_size-1-m2, result_real_neg, result_imag_neg);
				}
			}
			
			Two_point_objects.push_back(current_two_point);
		}

		current_time = std::time(nullptr);
		std::cout << current_time-start << "s: " << "Two-point objects have been computed.\n"; 

		// Compute ensemble mean of Fourier-Bessel coefficients
		mean(Two_point_objects, TwoPointFunction, TwoPointFunction_err);


	}

	// // Method to compute a two-point correlation function of the form <e_00 eml> (e_00 being the background coefficient)
	// void getTwoPointFunction_background(size_type centrality_index, complex_matrix<number_type> & TwoPointFunction, complex_matrix<number_type> & TwoPointFunction_err, std::time_t start)
	// {
	// 	std::time_t current_time = std::time(nullptr);

	// 	/* Bessel-decompose profile parts with quantum number m */
	// 	size_type mMax = (TwoPointFunction.rowsize()-1)/2; // TwoPoint contains m=-mMax...0...mMax
	// 	size_type lMax = TwoPointFunction.colsize();



	// 	// compute necessary bessel zeros
	// 	get_Bessel_deriv_zeros(mMax+1, lMax);

	// 	// compute two-point object e_00 e_ml for each profile in the centrality class

	// 	std::vector<complex_matrix<number_type>> Two_point_objects(0);

	// 	for (size_type k = 0; k < profiles_.size(); ++k)
	// 	{
	// 		// Skip profile if not in centrality class
	// 		if (!is_in_centrality_class(k, centrality_index))
	// 		{
	// 			continue;
	// 		}


	// 		complex_matrix<number_type> current_coeff(mMax+1, lMax);
	// 		BesselDecomposeProfile(k, centrality_index, current_coeff);
	// 		number_type current_background = getBackgroundCoefficient(k);
	// 		// fill current two-point object e_ml1 e*_ml2
	// 		complex_matrix<number_type> current_two_point(2*mMax+1, lMax);
	// 		for (size_type m = 0; m <= mMax; ++m)
	// 		{
	// 			for (size_type l = 0; l < lMax; ++l)
	// 			{
	// 				// Case m >= 0
	// 				number_type eml_pos_real = current_coeff.get_real(m, l);
	// 				number_type eml_pos_imag = current_coeff.get_imag(m, l);
					
	// 				number_type result_pos_real = current_background*eml_pos_real;
	// 				number_type result_pos_imag = current_background*eml_pos_imag;

	// 				current_two_point.set_entry(m+mMax, l, result_pos_real, result_pos_imag);

	// 				// case m < 0
	// 				number_type eml_neg_real = current_coeff.get_real(m, l);
	// 				number_type eml_neg_imag = current_coeff.get_imag(m, l)*(-1.0);
	// 				if (m%2 == 1)
	// 				{
	// 					eml_neg_real *= -1.0;
	// 					eml_neg_imag *= -1.0;
	// 				}
					
	// 				number_type result_neg_real = current_background*eml_neg_real;
	// 				number_type result_neg_imag = current_background*eml_neg_imag;

	// 				current_two_point.set_entry(mMax-m, l, result_neg_real, result_neg_imag);

	// 			}
	// 		}
			
	// 		Two_point_objects.push_back(current_two_point);
	// 	}

	// 	current_time = std::time(nullptr);
	// 	std::cout << current_time-start << "s: " << "Background two-point objects have been computed.\n"; 

	// 	// Compute ensemble mean of Fourier-Bessel coefficients
	// 	mean(Two_point_objects, TwoPointFunction, TwoPointFunction_err);


	// }

	// Bessel-decompose profile k and save data to current_coeff matrix
	void BesselDecomposeProfile(size_type k, size_type centrality_index, complex_matrix<number_type> & current_coeff)
	{
		size_type mMax = current_coeff.rowsize()-1;
		size_type lMax = current_coeff.colsize();	

		for (size_type m = 0; m <= mMax; ++m)
		{
			for (size_type l = 1; l <= lMax; ++l)
			{
				number_type real = integ_fourier_r(k, centrality_index, m, l, 0, R_-grid_step_);
				// choose correct back transformation formula
				if (m == 0 && l == 1)
				{
					real /= 1.0;
				}
				else if (m == 0 && l > 1)
				{
					real /= c(m, l-1); 
				}
				else
				{
					real /= c(m, l);
				}
				number_type imag;
				if (m == 0)
				{
					imag = 0;
				}
				else
				{
					imag = integ_fourier_r(k, centrality_index, Nm_-m, l, 0, R_-grid_step_)/c(m, l);
				}
				current_coeff.set_entry(m, l-1, real, imag);
			}
		}
	}



	// Bessel-decompose profile k with respect to mode m and save data to current_coeff matrix of the form 1 x lMax
	void BesselDecomposeProfile_fixed_m(size_type k, size_type centrality_index, size_type m, complex_matrix<number_type> & current_coeff)
	{
		assert(current_coeff.rowsize() == 1);
		size_type lMax = current_coeff.colsize();	

		
		for (size_type l = 1; l <= lMax; ++l)
		{
			number_type real = integ_fourier_r(k, centrality_index, m, l, 0, R_-grid_step_);
			// choose correct back transformation formula
			if (m == 0 && l == 1)
			{
				real /= 1.0;
			}
			else if (m == 0 && l > 1)
			{
				real /= c(m, l-1); 
			}
			else
			{
				real /= c(m, l);
			}
			number_type imag;
			if (m == 0)
			{
				imag = 0;
			}
			else
			{
				imag = integ_fourier_r(k, centrality_index, Nm_-m, l, 0, R_-grid_step_)/c(m, l);
			}
			current_coeff.set_entry(0, l-1, real, imag);
		}

	
	}

	// Bessel-decompose profile k with respect to mode l and save data to current_coeff matrix of the form m_size x 1
	void BesselDecomposeProfile_fixed_l(size_type k, size_type centrality_index, size_type l, complex_matrix<number_type> & current_coeff)
	{
		assert(current_coeff.colsize() == 1);
		size_type m_size = current_coeff.rowsize();	

		
		for (size_type m = 0; m < m_size; ++m)
		{
			number_type real = integ_fourier_r(k, centrality_index, m, l, 0, R_-grid_step_);
			// choose correct back transformation formula
			if (m == 0 && l == 1)
			{
				real /= 1.0;
			}
			else if (m == 0 && l > 1)
			{
				real /= c(m, l-1); 
			}
			else
			{
				real /= c(m, l);
			}
			number_type imag;
			if (m == 0)
			{
				imag = 0;
			}
			else
			{
				imag = integ_fourier_r(k, centrality_index, Nm_-m, l, 0, R_-grid_step_)/c(m, l);
			}
			current_coeff.set_entry(m, 0, real, imag);
		}

	
	}

	// // get coefficient e_00 for profile k
	// number_type getBackgroundCoefficient(size_type k)
	// {
	// 	return -1. + integ_r(k, 0, R_-grid_step_);
	// }

	/* Do a FFT decomposition of energy density profiles with respect to phi
	   yielding profiles of the form
	   E(m, r)
	*/
	void decompose_azimuthal()
	{
		for (size_type k = 0; k < profiles_.size(); ++k) // loop over profiles
		{
			// Set up matrix where (m, r)-profile will be saved
			gsl_matrix* profile = gsl_matrix_alloc(Nm_, Nr_);

			for (size_type i = 0; i < Nr_; ++i) // loop over radii
			{
				number_type radius = R_ * i / Nr_;
				number_type* e_m_r;
				e_m_r = new number_type[Nm_];
				// fill e_m_r with Energy(phi, r) values
				for (size_type j = 0; j < Nm_; ++j)
				{
					number_type phi = 2.*3.1415926 * j / Nm_;
					e_m_r[j] = interpolate(k, to_x(radius, phi), to_y(radius, phi));
				}
				// FFT it
				gsl_fft_real_radix2_transform(e_m_r, 1, Nm_);

				// save result in matrix
				for (size_type j = 0; j < Nm_; ++j)
				{
					gsl_matrix_set(profile, j, i, e_m_r[j]/Nm_);
				}
			}

			// save computed profile internally
			m_profiles_.push_back(profile);

			// free xy-interpolation objects
			gsl_spline2d_free(splines_[k]);
			gsl_interp_accel_free(xacc_[k]);
			gsl_interp_accel_free(yacc_[k]);

		}

		/* 
		Compute ensemble average of zero-th mode with yields the W(r) weight function needed for the
		Bessel decomposition -- for each centrality class
		*/

		for (size_type c = 1; c < percentiles_.size(); ++c)
		{
			gsl_vector* radii = gsl_vector_alloc(Nr_);
			gsl_vector* e_mean = gsl_vector_alloc(Nr_);
			gsl_vector* e_err = gsl_vector_alloc(Nr_);
			average_fourier(0, c, radii, e_mean, e_err);

			number_type* W_sites;
			W_sites = new number_type[Nr_];
			for (size_type i = 0; i < Nr_; ++i)
			{
				W_sites[i] = gsl_vector_get(e_mean, i)*3.1415926;
			}

			gsl_vector_free(radii);
			gsl_vector_free(e_mean);
			gsl_vector_free(e_err);

			// Initialize interpolation objects for W(r) interpolation
			const gsl_interp_type* interpolator = gsl_interp_cspline;
			W_interpolators_.push_back(interpolator);
			gsl_spline* spline = gsl_spline_alloc(W_interpolators_[c-1], Nr_);
			W_splines_.push_back(spline);
			gsl_interp_accel* acc = gsl_interp_accel_alloc();
			W_acc_.push_back(acc);
			gsl_spline_init(W_splines_[c-1], r_sites_, W_sites, Nr_);
		}

		// /*
		//    Subtract background function from the Fourier-decomposed data
		// */

		// for (size_type k = 0; k < m_profiles_.size(); ++k)
		// {
		// 	// Find out to which centrality class this profile belongs
		// 	for (size_type c = 1; c < percent_mult_.size(); ++c)
		// 	{
		// 		if(is_in_centrality_class(k, c))
		// 		{
		// 			// subtract correct background function from m=0 mode (actually not necessary since adding multiples of W(r) leaves coefficients invariant)
		// 			for (size_type l = 0; l < m_profiles_[k]->size2; ++l)
		// 			{
		// 				number_type current_value = gsl_matrix_get(m_profiles_[k], 0, l);
		// 				number_type radius = 1.0*(R_) * l / (Nr_+1);
		// 				gsl_matrix_set(m_profiles_[k], 0, l, current_value - 0.0*W(radius, c)/3.1415926 );
		// 			}
		// 			break;
		// 		}
		// 		else
		// 		{
		// 			continue;
		// 		}
		// 	}

		// }

		// // Sanity check: compute <e_0> again to see if it vanishes


		// for (size_type c = 1; c < percentiles_.size(); ++c)
		// {
		// 	gsl_vector* radii = gsl_vector_alloc(Nr_);
		// 	gsl_vector* e_mean = gsl_vector_alloc(Nr_);
		// 	gsl_vector* e_err = gsl_vector_alloc(Nr_);
		// 	average_fourier(0, c, radii, e_mean, e_err);
		// 	for (size_type i = 0; (i < e_mean->size) && (c == 1); ++i)
		// 		std::cout << gsl_vector_get(e_mean, i) << "\n";
		// }
		

	}



	

	// generate output file with values of W(r) for each centrality class
	void print_W(std::string filename, size_type Nr)
	{
		std::vector<std::vector<number_type>> Ws(percentiles_.size());
		for (size_type i = 0; i < Nr; ++i) // loop over radii
		{
			number_type radius = (R_-grid_step_) * i / Nr;
			(Ws[0]).push_back(radius);
			for (size_type c = 1; c < Ws.size(); ++c) // loop over centrality classes
			{
				(Ws[c]).push_back(W(radius, c));
			}
		}

		to_file(filename, Ws);
	}

	// generate output file with values of W(r) for each centrality class
	void print_W(std::string filename_begin)
	{
		// Ws structure: radius, W (sorted by central. class)
		std::vector<gsl_vector*> Ws(percentiles_.size());
		// dWs structure: radius, dW (sorted by central. class)
		std::vector<gsl_vector*> dWs(percentiles_.size());

		for (size_type c = 1; c < percentiles_.size(); ++c) // loop over centrality classes
		{
			// Compute <e_0>
			gsl_vector* radii = gsl_vector_alloc(Nr_);
			gsl_vector* e_mean = gsl_vector_alloc(Nr_);
			gsl_vector* e_err = gsl_vector_alloc(Nr_);

			average_fourier(0, c, radii, e_mean, e_err);


			// Scale them to get W(r)
			gsl_vector_scale(e_mean, 3.1415926);
			gsl_vector_scale(e_err, 3.1415926);


			// Save the at the correct position of Ws
			if (c == 1)
			{
				Ws[0] = radii;
				dWs[0] = radii;
			}


			Ws[c] = e_mean;
			dWs[c] = e_err;

		}

		// Save files

		to_file(filename_begin + ".txt" , Ws);
		to_file(filename_begin + "_error" + ".txt" , dWs);
		
	}

	// generate output file with <epsilon_m(r)> and uncertainty
	void print_e_m(std::string filename, size_type m, size_type c)
	{
		std::vector<gsl_vector*> e_ms(3);
	
		gsl_vector* radii = gsl_vector_alloc(Nr_);
		gsl_vector* e_mean = gsl_vector_alloc(Nr_);
		gsl_vector* e_err = gsl_vector_alloc(Nr_);
		average_fourier(m, c, radii, e_mean, e_err);

		e_ms[0] = radii;
		e_ms[1] = e_mean;
		e_ms[2] = e_err;

		to_file(filename, e_ms);

		gsl_vector_free(radii);
		gsl_vector_free(e_mean);
		gsl_vector_free(e_err);
		
	}

	

	// generate output file with values of W(r) for each centrality class
	void print_rho(std::string filename, size_type Nr)
	{
		std::vector<std::vector<number_type>> rhos(percentiles_.size());
		for (size_type i = 0; i < Nr; ++i) // loop over radii
		{
			number_type radius = (R_-grid_step_) * i / Nr;
			(rhos[0]).push_back(radius);
			for (size_type c = 1; c < rhos.size(); ++c) // loop over centrality classes
			{
				(rhos[c]).push_back(rho(radius, c));
			}
		}

		to_file(filename, rhos);
	}

	bool is_in_centrality_class(size_type profile_index, size_type centrality_index) const
	{
		assert(centrality_index != 0);
		number_type mult_min = percent_mult_[centrality_index];
		number_type mult_max = percent_mult_[centrality_index-1];
		number_type mult = multiplicities_[profile_index];
		if (mult <= mult_min || mult > mult_max ) // if profile not in centrality class
		{
			return false;
		}
		else 
		{
			return true;
		}
	}

	// compute average of Fourier decomposed data, specifically <E_m(r)> for given centrality class
	void average_fourier(size_type m, size_type centrality_index, gsl_vector* radii, gsl_vector* e_mean, gsl_vector* e_err)
	{
		assert(e_mean->size == Nr_);
		for (size_type i = 0; i < Nr_; ++i) // loop over radii
		{
			number_type r = R_ * i / Nr_; // current radius
			gsl_vector_set(radii, i, r);
			std::vector<number_type> e_ensemble(0);
			number_type mean_val;
			number_type mean_err;
			for (size_type k = 0; k < m_profiles_.size(); ++k) // loop over profiles
			{
				// Skip profile if not in centrality class
				if (!is_in_centrality_class(k, centrality_index))
				{
					continue;
				}

				number_type real_part = gsl_matrix_get(m_profiles_[k], m, i);
				number_type imag_part;
				if (m == 0)
				{
					imag_part = 0;
				}
				else
				{
					imag_part = gsl_matrix_get(m_profiles_[k], Nm_ - m, i);
				}

				number_type module = sqrt(real_part*real_part+imag_part*imag_part)*0 + real_part;

				
				//module *= 1.0*sign(real_part); // I want the actual value for W(r)
				

				e_ensemble.push_back(module);
			}
			mean(e_ensemble, mean_val, mean_err);
			gsl_vector_set(e_mean, i, mean_val);
			gsl_vector_set(e_err, i, mean_err);

		}	
	}

	// print data matrix k
	void print(size_type k) const
	{
		assert(k < profiles_.size());
		for (size_type i = 0; i < profiles_[k]->size1; ++i)
		{
			for (size_type j = 0; j < profiles_[k]->size2; ++j)
			{
				std::cout << gsl_matrix_get(profiles_[k], i, j) << " ";
			}
			std::cout << "\n";
		}
	}




private:
	// grid parameters
	number_type grid_max_;
	number_type grid_step_;
	size_type grid_size_;
	// vector of energy density profiles
	std::vector<gsl_matrix*> profiles_;

	// Random reaction plane
	bool random_reaction_plane_;
	std::vector<number_type> angles_;


	// vector of Fourier-decomposed (m, r)-profiles
	std::vector<gsl_matrix*> m_profiles_;

	// vector of impact_parameters
	std::vector<number_type> impact_parameters_;
	//vector of wounded nucleon number
	std::vector<number_type> n_participants_;
	//vector of multiplicities
	std::vector<number_type> multiplicities_;

	// vector of percentiles of interest
	std::vector<number_type> percentiles_;
	// corresponding multiplicities
	std::vector<number_type> percent_mult_;

	// x-y-interpolation objects
	number_type* x_sites_;
	number_type* y_sites_;
	std::vector<number_type*> z_profiles_;
	std::vector<const gsl_interp2d_type*> interpolators_;
	std::vector<gsl_spline2d*> splines_;
	std::vector<gsl_interp_accel*> xacc_;
	std::vector<gsl_interp_accel*> yacc_;

	// r-interpolation objects
	number_type* r_sites_;
	number_type R_;
	const size_type Nr_; // number of lattice sites in r direction
	const size_type Nm_; // number of m-modes for Fourier decomposition
	std::vector<const gsl_interp_type*> r_interpolators_;
	std::vector<gsl_spline*> r_splines_;
	std::vector<gsl_interp_accel*> racc_;

	// W(r) interpolation
	std::vector<const gsl_interp_type*> W_interpolators_;
	std::vector<gsl_spline*> W_splines_;
	std::vector<gsl_interp_accel*> W_acc_;
	std::vector<number_type> normalizations_;

	// Miscellaneous
	gsl_matrix* bessel_deriv_zeros_;




};

// struct to feed profile index and radius to phi integration function
template<typename size_type, typename number_type>
struct interpolate_at_r_params
{
	size_type profile_index;
	number_type radius;
	Collision<number_type>* pt_Collision;
};

// gsl function to call to integrate over the energy density with respect to phi at fixed R
double interpolate_at_r_wrapper(double phi, void* params)
{
	interpolate_at_r_params<std::size_t, double>* params_ = (interpolate_at_r_params<std::size_t, double>*) params;
	std::size_t k = params_->profile_index;
	double r = params_->radius;
	double x = to_x(r, phi);
	double y = to_y(r, phi);
	return params_->pt_Collision->interpolate(k, x, y);
}


// integrate over phi at fixed r with gsl integration routine
// Collision method, defined after (!) interpolate_at_r_params struct
template<class number_type>
number_type Collision<number_type>::integ_phi(std::size_t k, number_type r, number_type phi_lower, number_type phi_upper)
{

	//std::cout << "integ_phi\n";

	typedef std::size_t size_type;
	gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);
	number_type integral;
	number_type error;
	interpolate_at_r_params<size_type, number_type> params = {k, r, this};

	gsl_function F;
	F.function = &interpolate_at_r_wrapper;
	F.params = &params;
	gsl_integration_qags(&F, phi_lower, phi_upper, 1e-5, 1e-5, 1000, w, &integral, &error);
	//std::cout << " Radius: " << r << " Index: " << k <<  " Integral: " << integral << " Error abs: " << error << " Error rel: " << error/integral << "\n"; 
	gsl_integration_workspace_free(w);
	return integral;
}

// struct to feed profile index and mode number m to r integration function
template<typename size_type, typename number_type>
struct interpolate_at_m_params
{
	size_type profile_index;
	size_type mode_m;
	size_type mode_l;
	size_type centrality_index;
	Collision<number_type>* pt_Collision;
};

// gsl function to call to integrate over the mth Fourier mode of the energy density with respect to r at fixed m
double interpolate_at_m_wrapper(double r, void* params)
{
	interpolate_at_m_params<std::size_t, double>* params_ = (interpolate_at_m_params<std::size_t, double>*) params;
	std::size_t k = params_->profile_index;
	std::size_t m = params_->mode_m;
	std::size_t l = params_->mode_l;
	std::size_t c = params_->centrality_index;
	size_type Nm = params_->pt_Collision->Nm();

	// if (m < Nm/2)
	// {
	// 	return ((params_->pt_Collision->interpolate_fourier(k, m, r)) )*r*(params_->pt_Collision->Bessel(m, l, r, c)); // times r because of dr measure
	// }
	// else // m> Nm/2 calls refer to the imaginary part of E_(Nm-m) (r). 
	// {
	// 	return ((params_->pt_Collision->interpolate_fourier(k, m, r))  )*r*(params_->pt_Collision->Bessel(Nm-m, l, r, c)); // times r because of dr measure
	// }

	if (m < Nm/2)
	{
		if (m == 0 && l == 1)
		{
			return 2.*((params_->pt_Collision->interpolate_fourier(k, m, r)) )*r; // times r because of dr measure	
		}
		else if (m == 0 && l > 1)
		{
			return ((params_->pt_Collision->interpolate_fourier(k, m, r)) )*r*(params_->pt_Collision->Bessel(m, l-1, r, c)); // times r because of dr measure	
		}
		else
		{
			return ((params_->pt_Collision->interpolate_fourier(k, m, r)) )*r*(params_->pt_Collision->Bessel(m, l, r, c)); // times r because of dr measure	
		}

		
	}
	else // m> Nm/2 calls refer to the imaginary part of E_(Nm-m) (r). No case distinction needed because imaginary part = 0 for m=0
	{
		return ((params_->pt_Collision->interpolate_fourier(k, m, r))  )*r*(params_->pt_Collision->Bessel(Nm-m, l, r, c)); // times r because of dr measure
	}




	// if (m < Nm/2)
	// {
	// 	if (m == 0 && l == 1)
	// 	{
	// 		return (params_->pt_Collision->interpolate_fourier(k, m, r))*r; // times r because of dr measure
	// 	}
	// 	else if (m == 0 && l > 1)
	// 	{
	// 		return (params_->pt_Collision->interpolate_fourier(k, m, r))*r*(params_->pt_Collision->Bessel(m, l-1, r, c));
	// 	}
	// 	else
	// 	{
	// 		return (params_->pt_Collision->interpolate_fourier(k, m, r))*r*(params_->pt_Collision->Bessel(m, l, r, c)); // times r because of dr measure
	// 	}
	// }
	// else // m> Nm/2 calls refer to the imaginary part of E_(Nm-m) (r). No case distinction needed because imaginary part = 0 for m=0
	// {
	// 	return (params_->pt_Collision->interpolate_fourier(k, m, r))*r*(params_->pt_Collision->Bessel(Nm-m, l, r, c)); // times r because of dr measure
	// }	
	
}

// integrate over r at fixed Fourier-mode m and Bessel(m,l) as weight with gsl integration routine
// Collision method, defined after (!) interpolate_at_m_params struct
template<class number_type>
number_type Collision<number_type>::integ_fourier_r(size_type k, size_type centrality_index, size_type m, size_type l, number_type r_lower, number_type r_upper)
{

	//std::cout << "integ_fourier_r\n";

	typedef std::size_t size_type;
	gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);
	number_type integral;
	number_type error;
	interpolate_at_m_params<size_type, number_type> params = {k, m, l, centrality_index, this};

	gsl_function F;
	F.function = &interpolate_at_m_wrapper;
	F.params = &params;
	gsl_integration_qags(&F, r_lower, r_upper, 1e-5, 1e-5, 1000, w, &integral, &error);
	//std::cout << " Radius: " << r << " Index: " << k <<  " Integral: " << integral << " Error abs: " << error << " Error rel: " << error/integral << "\n"; 
	gsl_integration_workspace_free(w);
	return integral;
}


// integrate dr r  dphi epsilon(r, phi)

// gsl function to call to integrate over the mth Fourier mode of the energy density with respect to r at fixed m
double interpolate_r_wrapper(double r, void* params)
{
	interpolate_at_m_params<std::size_t, double>* params_ = (interpolate_at_m_params<std::size_t, double>*) params;
	std::size_t k = params_->profile_index;

	return ((params_->pt_Collision->interpolate_fourier(k, 0, r)) )*r*2.*3.1415926; 
}
template<class number_type>
number_type Collision<number_type>::integ_r(size_type k, number_type r_lower, number_type r_upper)
{

	//std::cout << "integ_r\n";

	typedef std::size_t size_type;
	gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);
	number_type integral;
	number_type error;
	interpolate_at_m_params<size_type, number_type> params = {k, 0, 1, 1, this};

	gsl_function F;
	F.function = &interpolate_r_wrapper;
	F.params = &params;
	gsl_integration_qags(&F, r_lower, r_upper, 1e-5, 1e-5, 1000, w, &integral, &error);
	//std::cout << " Integral: " << integral << " Error abs: " << error << " Error rel: " << error/integral << "\n"; 
	gsl_integration_workspace_free(w);
	return integral;
}

// integrate dr r  W(r)

// gsl function to call to integrate over the mth Fourier mode of the energy density with respect to r at fixed m
double interpolate_r_Wr_wrapper(double r, void* params)
{
	interpolate_at_m_params<std::size_t, double>* params_ = (interpolate_at_m_params<std::size_t, double>*) params;
	std::size_t c = params_->centrality_index;

	return (params_->pt_Collision->W(r, c)) *r; 
}
template<class number_type>
number_type Collision<number_type>::integ_r_Wr(size_type centrality_index, number_type r_lower, number_type r_upper)
{

	//std::cout << "integ_r_Wr\n";

	typedef std::size_t size_type;
	gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);
	number_type integral;
	number_type error;
	interpolate_at_m_params<size_type, number_type> params = {0, 0, 1, centrality_index, this};

	gsl_function F;
	F.function = &interpolate_r_Wr_wrapper;
	F.params = &params;
	gsl_integration_qags(&F, r_lower, r_upper, 1e-5, 1e-5, 1000, w, &integral, &error);
	//std::cout << " Radius: " << r << " Index: " << k <<  " Integral: " << integral << " Error abs: " << error << " Error rel: " << error/integral << "\n"; 
	gsl_integration_workspace_free(w);
	return integral;
}

// struct to feed Collision object to rho(r) function
template<typename size_type, typename number_type>
struct interpolate_W_params
{
	Collision<number_type>* pt_Collision;
	size_type centrality_index;
};

// gsl function to call for rho(r) function
double interpolate_W_wrapper(double r, void* params)
{
	interpolate_W_params<std::size_t, double>* params_ = (interpolate_W_params<std::size_t, double>*) params;
	return (params_->pt_Collision->W(r, params_->centrality_index))*r; // times r because of dr measure
}

// rho(r) function, see description above declaration in class
// Collision method, defined after (!) interpolate_W_params struct
template<class number_type>
number_type Collision<number_type>::rho(number_type r, size_type centrality_index)
{

	//std::cout << "integ_rho\n";

	typedef std::size_t size_type;
	gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);
	number_type integral;
	number_type error;
	interpolate_W_params<size_type, number_type> params = {this, centrality_index};

	gsl_function F;
	F.function = &interpolate_W_wrapper;
	F.params = &params;
	gsl_integration_qags(&F, 0, r, 1e-5, 1e-5, 1000, w, &integral, &error);
	gsl_integration_workspace_free(w);
	return sqrt(integral*2.0);
}


// Test function: integrate W(r)r J_0(z_l^0 rho(r))

// struct to feed Collision object to test function
template<typename size_type, typename number_type>
struct test_params
{
	Collision<number_type>* pt_Collision;
	size_type centrality_index;
};
// gsl function to call to integrate over the mth Fourier mode of the energy density with respect to r at fixed m
double test_wrapper(double r, void* params)
{
	test_params<std::size_t, double>* params_ = (test_params<std::size_t, double>*) params;
	std::size_t m = 1;
	std::size_t l = 2;
	std::size_t c = params_->centrality_index;
	params_->pt_Collision->get_Bessel_deriv_zeros(m+1, l+1);

	return (params_->pt_Collision->Bessel(m, l, r, c)) * r * (params_->pt_Collision->W(r, c));	
	
}

template<class number_type>
number_type Collision<number_type>::integ_test()
{

	//std::cout << "integ_test\n";

	typedef std::size_t size_type;
	gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);
	number_type integral;
	number_type error;
	size_type centrality_index = 4;
	test_params<size_type, number_type> params = {this, centrality_index};

	gsl_function F;
	F.function = &test_wrapper;
	F.params = &params;
	gsl_integration_qags(&F, 0, R_-grid_step_, 1e-5, 1e-5, 1000, w, &integral, &error);
	//std::cout << " Radius: " << r << " Index: " << k <<  " Integral: " << integral << " Error abs: " << error << " Error rel: " << error/integral << "\n"; 
	gsl_integration_workspace_free(w);
	return integral;
}

#endif