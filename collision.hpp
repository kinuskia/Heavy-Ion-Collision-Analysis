#ifndef COLLISION_HPP
#define COLLISION_HPP

#include "read_data.hpp"

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
	, data_(0)
	{}


	// read in a data file
	void read_in(std::string filename)
	{
		Matrix<number_type> data(grid_size_, grid_size_);
		read_data(filename, data, 8);
		data_.push_back(data);
	}

	// normalize matrix data (integral = 1)
	void normalize(number_type energy)
	{
		for (size_type k = 0; k < data_.size(); ++k)
		{
			// compute total energy
			number_type sum = 0;
			for (size_type i = 0; i < data_[k].rowsize(); ++i)
			{
				for (size_type j = 0; j < data_[k].colsize(); ++j)
				{
					sum += data_[k](i, j);
				}
			}
			sum *= grid_step_*grid_step_;
			// rescale
			for (size_type i = 0; i < data_[k].rowsize(); ++i)
			{
				for (size_type j = 0; j < data_[k].colsize(); ++j)
				{
					data_[k](i, j) *= energy/sum;
				}
			}
		}
	}

	// centralize matrix data (barycenter at N/2,N/2)
	void centralize()
	{
		for (size_type k = 0; k < data_.size(); ++k)
		{
			// compute barycentre
			number_type x_bary = 0;
			number_type y_bary = 0;
			for (size_type i = 0; i < data_[k].rowsize(); ++i)
			{
				for (size_type j = 0; j < data_[k].colsize(); ++j)
				{
					number_type x = j;
					number_type y = i;
					x_bary += data_[k](i, j)*x;
					y_bary += data_[k](i, j)*y;
				}
			}
			x_bary *= grid_step_*grid_step_;
			y_bary *= grid_step_*grid_step_;


			// make temporary copy
			Matrix<number_type> m_copy = data_[k];

			// center the barycentre
			size_type xB = round(x_bary);
			size_type yB = round(y_bary);
			int dx = m_copy.colsize()/2-xB;
			int dy = m_copy.rowsize()/2-yB;
			std::cout << "x: " << xB << " y: " << yB << "\n";
			
			for (int i = 0; i < m_copy.rowsize(); ++i)
			{
				for (int j = 0; j < m_copy.colsize(); ++j)
				{
					// deal with out-of bound indices
					if (i-dy < 0 || j-dx < 0 || i-dy >= m_copy.colsize() || j-dx >= m_copy.rowsize())
					{
						data_[k](i, j) = 0;
					}
					else
					{
						data_[k](i, j) = m_copy(i-dy, j-dx);
					}
				}
			}
		}
	}

	// print data matrix k
	void print(size_type k) const
	{
		assert(k < data_.size());
		for (size_type i = 0; i < data_[k].rowsize(); ++i)
		{
			for (size_type j = 0; j < data_[k].colsize(); ++j)
			{
				std::cout << data_[k](i, j) << " ";
			}
			std::cout << "\n";
		}
	}


private:
	number_type grid_max_;
	number_type grid_step_;
	size_type grid_size_;
	std::vector<Matrix<number_type>> data_;
};

#endif