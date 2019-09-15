#ifndef MODEL_HPP
#define MODEL_HPP

#include <iostream>
#include <cmath>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_spline2d.h>
#include "../auxiliary/read_data.hpp"

#include <fstream>

template<class REAL>
class Model
{
public:
	typedef std::size_t size_type;
	typedef REAL number_type;

	// constructor
	Model(number_type m, number_type Q0)
	: m_(m)
	, Q0_(Q0)
	, pi_(3.1415926)
	, gridmax_(10)
	, gridstep_(0.2)
	{}

	// initialize One-Point function
	void initialize_OnePoint(std::string filename, size_type N, const gsl_interp2d_type* interp_method, number_type gridmax, number_type gridstep)
	{
		gsl_matrix* profile = gsl_matrix_alloc(N, N);
		number_type dummy;
		read_data(filename, profile, dummy, dummy, dummy);
		

		const gsl_interp2d_type* interpolator = interp_method;
		OnePoint_interpolator_ = interpolator;
		number_type* za = new number_type[N*N];
		gsl_spline2d* spline = gsl_spline2d_alloc(interpolator, N, N);
		OnePoint_spline_ = spline;
		gsl_interp_accel* xacc = gsl_interp_accel_alloc();
		gsl_interp_accel* yacc = gsl_interp_accel_alloc();
		xacc_ = xacc;
		yacc_ = yacc;

		// set interpolation grid values
		number_type* x_sites = new number_type[N];
		number_type* y_sites = new number_type[N];
		gridmax_ = gridmax;
		gridstep_ = gridstep;
		number_type min = -gridmax_ +gridstep_*0.51 ;
		number_type max = gridmax_ -gridstep_*0.51 ;
		for (size_type i = 0; i < N ; ++i)
		{
			x_sites[i] = min + (max-min)*i/(N-1);
			y_sites[i] = x_sites[i];
		}

		for (size_type i = 0; i < N; ++i)
		{
			for (size_type j = 0; j < N; ++j)
			{
				gsl_spline2d_set(OnePoint_spline_, za, i, j, gsl_matrix_get(profile, N-1 - j, i));
			}
		}

		gsl_spline2d_init(OnePoint_spline_, x_sites, y_sites, za, N, N);

		// free profile
		gsl_matrix_free(profile);
	}




	// Define one-point position space correlation function
	// units in fm
	number_type OnePoint(number_type x, number_type y)
	{
		// number_type r = sqrt(x*x+y*y);
		// return W(r)/pi_;
		number_type result = gsl_spline2d_eval(OnePoint_spline_, x, y, xacc_, yacc_);
		if (result < 0)
		{
			result = -result;
		}

		return result;

	}


	// Access current event that has been read in
	number_type Event(Collision<number_type> Trento, size_type index, number_type x, number_type y)
	{	
		
		number_type result = Trento.interpolate(index, x, y);
		
		if (result < 0) //easy presciption if due to rounding errors the mean value is e.g. -1.3e-15
		{
			result = -result;
		}

		return result;
	}

	// test method to print One-Point profile
	void print_OnePoint(size_type N)
	{
		number_type min = -gridmax_ +gridstep_*0.51 ;
		number_type max = gridmax_ -gridstep_*0.51 ;
		
		for (size_type i = 0; i < N; ++i)
		{
			for (size_type j = 0; j < N; ++j)
			{
				number_type x = min + (max-min)*j/(N-1);
				number_type y = max - (max-min)*i/(N-1);
				number_type current = OnePoint(x, y);
				std::cout << current << " "; 
			}
			std::cout << "\n";
		}
	}

	// Define connected (!) position space two-point correlation function TwoPoint(x1, y1) = TwoPoint(x1, y1)*delta(x2-x1)*delta(y2-y1)
	number_type TwoPoint(number_type x, number_type y)
	{
		// number_type r = sqrt(x*x+y*y);
		// number_type alpha_s = 0.4095;
		// number_type g = sqrt(4.*pi_*alpha_s);
		// number_type W0 = W(0.);
		// number_type Q0 = 1.1*1+2.4*0; //GeV
		// number_type scale = Q0*Q0*Q0*Q0/3/alpha_s/W0;
		// number_type hc = 0.197327;
		// number_type Q2 = sqrt(OnePoint(x, y))*sqrt(3./4)*g*sqrt(scale);
		// return 32.*pi_/9/g/g/g/g*Q2*Q2*Q2*log(Q2/m_/m_)/scale/scale*hc*hc;


		//number_type r = sqrt(x*x+y*y);
		//number_type alpha_s = 0.4095;
		//number_type g = sqrt(4.*pi_*alpha_s);
		number_type p0 = OnePoint(0,0);
		number_type Q0 = Q0_; //GeV
		//number_type scale = Q0*Q0*Q0*Q0/3/alpha_s/W0;
		number_type hc = 0.197327;
		//number_type Q2 = sqrt(OnePoint(x, y))*sqrt(3./4)*g*sqrt(scale);

		number_type R_cutoff = 0.010;
		number_type r = sqrt(x*x+y*y);
		number_type phi = atan2(y,x);
        if (r < R_cutoff)
        {
        	x = R_cutoff*cos(phi);
        	y = R_cutoff*sin(phi);	
        }

		//return 2.*pi_*hc*hc*p0*p0/Q0/Q0*pow((OnePoint(x, y)/p0),3./2)*log(Q0*Q0/m_/m_*sqrt(OnePoint(x, y)/p0));
		return 2.*pi_*hc*hc*p0*p0/Q0/Q0*pow((Event(x, y)/p0),3./2)*log(Q0*Q0/m_/m_*sqrt(Event(x, y)/p0));
		
	}

private:

	const gsl_interp2d_type* OnePoint_interpolator_;
	gsl_spline2d* OnePoint_spline_;
	gsl_interp_accel* xacc_;
	gsl_interp_accel* yacc_;
	number_type gridmax_;
	number_type gridstep_;
	

	//
	number_type m_;
	number_type Q0_;
	number_type pi_;
};


#endif