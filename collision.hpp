#ifndef COLLISION_HPP
#define COLLISION_HPP

#include "read_data.hpp"
#include "polar.hpp"
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


	// read in a data file
	void read_in(std::string filename)
	{
		Matrix<number_type> data(grid_size_, grid_size_);
		read_data(filename, data, 8);
		profiles_.push_back(data);
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


	// average over azimuthal angle, result (r_i, <E_i>, error)
	void average_azimuthal()
	{
		for (size_type k = 0; k < profiles_.size(); ++k)
		{
			//
		}
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