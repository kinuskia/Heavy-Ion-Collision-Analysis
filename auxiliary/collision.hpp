#ifndef COLLISION_HPP
#define COLLISION_HPP

#include "read_data.hpp"
#include "/Users/Kianusch/Documents/Numerical_analysis/auxiliary_files/storage.hpp"
#include "/Users/Kianusch/Documents/Numerical_analysis/auxiliary_files/vector.hpp"
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
	{}

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
		Matrix<number_type> data(grid_size_, grid_size_);
		read_data(filename, data, 8);
		profiles_.push_back(data);
	}

	// read in data files
	void read_in(std::string filename_begin, std::string fileformat, size_type n_files)
	{
		for (size_type k = 0; k < n_files; ++k)
		{
			std::string filename = get_filename(filename_begin, fileformat, k, n_files);
			Matrix<number_type> data(grid_size_, grid_size_);
			read_data(filename, data, 8);
			profiles_.push_back(data);
		}
	}

	// normalize matrix data (integral = 1)
	void normalize(number_type energy)
	{
		for (size_type k = 0; k < profiles_.size(); ++k)
		{
			// compute total energy
			number_type sum = 0;
			for (size_type i = 0; i < profiles_[k].rowsize(); ++i)
			{
				for (size_type j = 0; j < profiles_[k].colsize(); ++j)
				{
					sum += profiles_[k](i, j);
				}
			}
			sum *= grid_step_*grid_step_;
			// rescale
			for (size_type i = 0; i < profiles_[k].rowsize(); ++i)
			{
				for (size_type j = 0; j < profiles_[k].colsize(); ++j)
				{
					profiles_[k](i, j) *= energy/sum;
				}
			}
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
			for (size_type i = 0; i < profiles_[k].rowsize(); ++i)
			{
				for (size_type j = 0; j < profiles_[k].colsize(); ++j)
				{
					number_type x = j;
					number_type y = i;
					x_bary += profiles_[k](i, j)*x;
					y_bary += profiles_[k](i, j)*y;
				}
			}
			x_bary *= grid_step_*grid_step_;
			y_bary *= grid_step_*grid_step_;


			// make temporary copy
			Matrix<number_type> m_copy = profiles_[k];

			// center the barycentre
			int dx = round((1.0*m_copy.colsize()-1)/2-x_bary);
			int dy = round((1.0*m_copy.rowsize()-1)/2-y_bary);
			//std::cout << "x: " << xB << " y: " << yB << "\n";
			
			for (int i = 0; i < m_copy.rowsize(); ++i)
			{
				for (int j = 0; j < m_copy.colsize(); ++j)
				{
					// deal with out-of bound indices
					if (i-dy < 0 || j-dx < 0 || i-dy >= m_copy.colsize() || j-dx >= m_copy.rowsize())
					{
						profiles_[k](i, j) = 0;
					}
					else
					{
						profiles_[k](i, j) = m_copy(i-dy, j-dx);
					}
				}
			}
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

	// interpolate energy density profil i at position (x,y)
	// now: bilinear interpolation
	number_type interpolate(size_type k, number_type x, number_type y) const
	{
		// Step 1: Find the grid square of (x,y), meaning the four tabulated
		// points that surround (x,y): (I,J) -> (I+1, J+1)
		// Let us enumerate the four points counter-clockwise from the lower left
		number_type e1;
		number_type e2;
		number_type e3;
		number_type e4;
		number_type xL;
		number_type xR;
		number_type yT;
		number_type yB;
		size_type I;
		size_type J;
		for (size_type i = 0; i < profiles_[k].rowsize(); ++i)
		{
			number_type Y = y_from_index(i);
			if (Y < y)
			{
				yB = Y;
				yT = Y+grid_step_;
				I = i-1;
				break;
			}
		}
		for (size_type j = 0; j < profiles_[k].colsize(); ++j)
		{
			number_type X = x_from_index(j);
			if (X > x)
			{
				xR = X;
				xL = X-grid_step_;
				J = j-1;
				break;
			}

		}
		e1 = profiles_[k](I+1, J);
		e2 = profiles_[k](I+1, J+1);
		e3 = profiles_[k](I, J+1);
		e4 = profiles_[k](I, J);

		// Step 2: Estimate energy density at (x,y)
		// relative coordinates (t,u) with respect to the grid square
		number_type t = (x-xL)/grid_step_;
		number_type u = (y-yB)/grid_step_;

		number_type energy = e1*(1.-t)*(1.-u)+e2*(t)*(1.-u)+e3*(t)*(u)+e4*(1.-t)*(u);

		return energy;
	}


	// average over azimuthal angle save to data storage
	void average_azimuthal(Vector<number_type> & radii, Vector<number_type> & mean_vec, Vector<number_type> & mean_error)
	{
		size_type N = mean_vec.size();
		Storage<number_type> profiles_averaged(1);
		for (size_type k = 0; k < profiles_.size(); ++k)
		{
			number_type R = grid_max_-1;
			for (size_type l = 0; l < N; ++l)
			{
				number_type r = number_type(l)/N*R;
				if (k == 0)
				{
					radii[l] = r;
				}
				number_type e_average = integ_phi(k, r, 0., 2.*3.1415926)/(2.*3.1415926);
				profiles_averaged.read_in(e_average);
			}
		}
		profiles_averaged.set_n_variables(radii.size());
		profiles_averaged.mean(mean_vec, mean_error, false);
	}

	// print data matrix k
	void print(size_type k) const
	{
		assert(k < profiles_.size());
		for (size_type i = 0; i < profiles_[k].rowsize(); ++i)
		{
			for (size_type j = 0; j < profiles_[k].colsize(); ++j)
			{
				std::cout << profiles_[k](i, j) << " ";
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
	std::vector<Matrix<number_type>> profiles_;




};

#endif