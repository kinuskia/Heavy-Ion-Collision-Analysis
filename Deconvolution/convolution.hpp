#ifndef CONVOLUTION_HPP
#define CONVOLUTION_HPP

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_spline2d.h>
#include "../auxiliary/read_data.hpp"
#include "../auxiliary/to_file.hpp"


template <class REAL>
class Convolution
{
public:
	typedef std::size_t size_type;
	typedef REAL number_type;

	Convolution(number_type grid_max, number_type grid_step)
	: grid_max_(grid_max)
	, grid_step_(grid_step)
	, grid_size_(size_type(2.*grid_max/grid_step))
	{}

	// Functions to convert index pair (i,j) to cartesian coordinates (x,y)
	number_type x_from_index(size_type j) const
	{
		return -grid_max_+grid_step_/2 + grid_step_*j;
	}

	number_type y_from_index(size_type i) const
	{
		return grid_max_-grid_step_/2 - grid_step_*i;
	}


	// import convolved profile
	void get_convolved(std::string filename)
	{
		gsl_matrix* convolved = gsl_matrix_alloc(grid_size_, grid_size_);
		number_type dummy;
		read_data(filename, convolved, dummy, dummy, dummy);

		const gsl_interp2d_type* interp_method = gsl_interp2d_bicubic;
		const gsl_interp2d_type* interpolator = interp_method;
		OnePoint_interpolator_ = interpolator;
		number_type* za = new number_type[grid_size_*grid_size_];
		gsl_spline2d* spline = gsl_spline2d_alloc(interpolator, grid_size_, grid_size_);
		OnePoint_spline_ = spline;
		gsl_interp_accel* xacc = gsl_interp_accel_alloc();
		gsl_interp_accel* yacc = gsl_interp_accel_alloc();
		xacc_ = xacc;
		yacc_ = yacc;	

		// set interpolation grid values
		number_type* x_sites = new number_type[grid_size_];
		number_type* y_sites = new number_type[grid_size_];
	
		number_type min = x_from_index(0) +grid_step_*0.01 ;
		number_type max = x_from_index(grid_size_-1) -grid_step_*0.01 ;
		size_type N = grid_size_;
		for (size_type i = 0; i < N ; ++i)
		{
			x_sites[i] = min + (max-min)*i/(N-1);
			y_sites[i] = x_sites[i];
		}

			for (size_type i = 0; i < N; ++i)
			{
				for (size_type j = 0; j < N; ++j)
				{
					gsl_spline2d_set(OnePoint_spline_, za, i, j, gsl_matrix_get(convolved, N-1 - j, i));
				}
			}

		gsl_spline2d_init(OnePoint_spline_, x_sites, y_sites, za, N, N);

		// free profile
		gsl_matrix_free(convolved);


	}

	number_type interpolate_convolved(number_type x, number_type y)
	{
		return gsl_spline2d_eval(OnePoint_spline_, x, y, xacc_, yacc_);
	}

	// test method to print One-Point profile
	void print_OnePoint(std::string outfile, size_type N)
	{
		number_type min = -grid_max_ +grid_step_ *0.51 ;
		number_type max = grid_max_ -grid_step_*0.51 ;
		gsl_matrix* output = gsl_matrix_alloc(N, N);
		for (size_type i = 0; i < N; ++i)
		{
			for (size_type j = 0; j < N; ++j)
			{
				number_type x = min + (max-min)*j/(N-1);
				number_type y = max - (max-min)*i/(N-1);
				number_type current = convolution(x, y);
				gsl_matrix_set(output, i, j, current);
			}
		}

		to_file(outfile, output);
		gsl_matrix_free(output);
	}

	// convolution kernel
	number_type kernel(number_type x, number_type y)
	{
		number_type sigma = 1.6;
		number_type pi = 3.1415926;

		return 1./2/pi/sigma/sigma*exp(-1./2/sigma/sigma*(x*x+y*y));
	}

	// convolve
	number_type convolution(number_type x, number_type y)
	{
		number_type min = -grid_max_ +grid_step_ *0.51 ;
		number_type max = grid_max_ -grid_step_*0.51 ;
		number_type integral_xy = 0; 
		size_type N = 50;
		for (size_type i = 0; i <= N; ++i)
		{
			number_type y1 = min + (max-min)*i/N;
			number_type integral_x = 0;
			for (size_type j = 0; j <= N; ++j)
			{
				number_type x1 = min + (max-min)*j/N;
				number_type f = interpolate_convolved(x1, y1)*kernel(x-x1, y-y1);
				if ((j==0)||(j==N))
				{
					integral_x += f/2;
				}
				else
				{
					integral_x += f;
				}
			}
			integral_x *= (max-min)/N;

			if ((i==0)||(i==N))
			{
				integral_xy += integral_x/2;
			}
			else
			{
				integral_xy += integral_x;
			}
		}
		integral_xy *= (max-min)/N;

		return integral_xy;
	}


private:
	number_type grid_max_;
	number_type grid_step_;
	size_type grid_size_;

	const gsl_interp2d_type* OnePoint_interpolator_;
	gsl_spline2d* OnePoint_spline_;
	gsl_interp_accel* xacc_;
	gsl_interp_accel* yacc_;


};


#endif