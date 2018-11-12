#ifndef COLLISION_HPP
#define COLLISION_HPP

#include <vector>
#include "read_data.hpp"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>
#include <algorithm>

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
	, z_profiles_(0)
	, interpolators_(0)
	, splines_(0)
	, xacc_(0)
	, yacc_(0) 
	{
		// initialize grid sites
		x_sites_ = new number_type[grid_size_];
		y_sites_ = new number_type[grid_size_];
		for (size_type i = 0; i < grid_size_; ++i)
		{
			x_sites_[i] = x_from_index(i);
			y_sites_[i] = -y_from_index(i);
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

	// read in single file
	void read_in(std::string filename)
	{
		gsl_matrix* data = gsl_matrix_alloc(grid_size_, grid_size_);
		read_data(filename, data, 8);
		profiles_.push_back(data);
	}

	// free profiles
	void free()
	{
		for (size_type i = 0; i < profiles_.size(); ++i)
		{
			gsl_matrix_free(profiles_[i]);
		}
	}

	// read in data files
	void read_in(std::string filename_begin, std::string fileformat, size_type n_files)
	{
		for (size_type k = 0; k < n_files; ++k)
		{
			std::string filename = get_filename(filename_begin, fileformat, k, n_files);
			gsl_matrix* data = gsl_matrix_alloc(grid_size_, grid_size_);
			read_data(filename, data, 8);
			profiles_.push_back(data);
		}
	}


			// initialize interpolation
	void initialize_interpolation(const gsl_interp2d_type* interp_method)
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
			
			// rescale
			gsl_matrix_scale(profiles_[k], energy/sum);
		}
	}

	// centralize matrix data (barycenter at N/2,N/2)
	void centralize()
	{
		for (size_type k = 0; k < profiles_.size(); ++k)
		{
			// compute barycentre
			number_type x_bary = 0;
			number_type y_bary = 0;
			for (size_type i = 0; i < profiles_[k]->size1; ++i)
			{
				for (size_type j = 0; j < profiles_[k]->size2; ++j)
				{
					number_type x = j;
					number_type y = i;
					x_bary += gsl_matrix_get(profiles_[k], i, j)*x;
					y_bary += gsl_matrix_get(profiles_[k], i, j)*y;
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

	// Functions to convert index pair (i,j) to cartesian coordinates (x,y)
	number_type x_from_index(size_type j) const
	{
		return -grid_max_+grid_step_/2 + grid_step_*j;
	}

	number_type y_from_index(size_type i) const
	{
		return grid_max_-grid_step_/2 - grid_step_*i;
	}

	// Functions to transform polar coordinates to cartesian coordinates
	number_type to_x(number_type r, number_type phi) const
	{
		return r*cos(phi);
	}
	number_type to_y(number_type r, number_type phi) const
	{
		return r*sin(phi);
	}

	// methods to integrate over phi (trapezoid rule) the profil k
	number_type integ_phi(size_type k, number_type r, number_type phi_lower, number_type phi_upper)
	{
		size_type N = grid_size_;
		number_type h = (phi_upper - phi_lower)/N;
		number_type integral = interpolate(k, to_x(r, phi_lower), to_y(r, phi_lower));
		integral += interpolate(k, to_x(r, phi_upper), to_y(r, phi_upper));
		integral *= 0.5;
		for (size_type l = 1; l < N; ++l)
		{
			number_type phi = phi_lower + h*l;
			integral += interpolate(k, to_x(r, phi), to_y(r, phi));
		}
		integral *= h;

		return integral;

	}

	// // interpolate energy density profil i at position (x,y)
	// // now: bilinear interpolation
	// number_type interpolate(size_type k, number_type x, number_type y) const
	// {
	// 	// Step 1: Find the grid square of (x,y), meaning the four tabulated
	// 	// points that surround (x,y): (I,J) -> (I+1, J+1)
	// 	// Let us enumerate the four points counter-clockwise from the lower left
	// 	number_type e1;
	// 	number_type e2;
	// 	number_type e3;
	// 	number_type e4;
	// 	number_type xL;
	// 	number_type xR;
	// 	number_type yT;
	// 	number_type yB;
	// 	size_type I;
	// 	size_type J;
	// 	for (size_type i = 0; i < profiles_[k]->size1; ++i)
	// 	{
	// 		number_type Y = y_from_index(i);
	// 		if (Y < y)
	// 		{
	// 			yB = Y;
	// 			yT = Y+grid_step_;
	// 			I = i-1;
	// 			break;
	// 		}
	// 	}
	// 	for (size_type j = 0; j < profiles_[k]->size2; ++j)
	// 	{
	// 		number_type X = x_from_index(j);
	// 		if (X > x)
	// 		{
	// 			xR = X;
	// 			xL = X-grid_step_;
	// 			J = j-1;
	// 			break;
	// 		}

	// 	}
	// 	e1 = gsl_matrix_get(profiles_[k], I+1, J);
	// 	e2 = gsl_matrix_get(profiles_[k], I+1, J+1);
	// 	e3 = gsl_matrix_get(profiles_[k], I, J+1);
	// 	e4 = gsl_matrix_get(profiles_[k], I, J);

	// 	// Step 2: Estimate energy density at (x,y)
	// 	// relative coordinates (t,u) with respect to the grid square
	// 	number_type t = (x-xL)/grid_step_;
	// 	number_type u = (y-yB)/grid_step_;

	// 	number_type energy = e1*(1.-t)*(1.-u)+e2*(t)*(1.-u)+e3*(t)*(u)+e4*(1.-t)*(u);

	// 	return energy;
	// }

	number_type interpolate(size_type k, number_type x, number_type y) const
	{
		return gsl_spline2d_eval(splines_[k], x, y, xacc_[k], yacc_[k]);
	}


	// average over azimuthal angle save to data storage
	void average_azimuthal(gsl_vector* radii, gsl_vector* mean_vec, gsl_vector* mean_error)
	{
		size_type N = radii->size;

		// For each profile compute e_average as function of r
		number_type R = grid_max_-1;
		for (size_type l = 0; l < N; ++l) // loop over radii
		{
			number_type r = number_type(l)/N*R;
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

	// interpolation objects
	number_type* x_sites_;
	number_type* y_sites_;
	std::vector<number_type*> z_profiles_;
	std::vector<const gsl_interp2d_type*> interpolators_;
	std::vector<gsl_spline2d*> splines_;
	std::vector<gsl_interp_accel*> xacc_;
	std::vector<gsl_interp_accel*> yacc_;


};

#endif