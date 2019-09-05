#ifndef MODEL_HPP
#define MODEL_HPP

#include <iostream>
#include <cmath>

#include <gsl/gsl_spline.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_spline2d.h>
#include <gsl/gsl_sf_bessel.h>

#include "../auxiliary/read_data.hpp"



#include <fstream>

template<class REAL>
class Model
{
public:
	typedef std::size_t size_type;
	typedef REAL number_type;

	// constructor
	Model(number_type m)
	: m_(m)
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

	// Define connected (!) position space two-point correlation function
	number_type TwoPoint(number_type x1, number_type y1, number_type x2, number_type y2)
	{
        number_type nc =3.;
        number_type hc = .1973;
        number_type r = sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
        r /= hc;

        number_type R_cutoff = 0.010/hc;
        if (r < R_cutoff)
        {
        	r = R_cutoff;	
        }

        number_type b1= (x1+x2)/2.;
        number_type b2= (y1+y2)/2.;
        number_type casimir= (nc*nc-1.)/2./nc;

        number_type alpha_s = 0.4095*0+0.25;
		number_type g = sqrt(4.*pi_*alpha_s);
		number_type p0 = OnePoint(0., 0.);
		number_type Q0 = 1.1*0+1.24; //GeV
		number_type scale = Q0*Q0*Q0*Q0/3./pi_/alpha_s/p0;
		number_type Q2 = Q0*Q0*sqrt(OnePoint(b1,b2)/p0);
       
        number_type coupling=g;
        number_type infrared_Regulator= m_*1e-2; //0-1 1e-2
        number_type modifiedGamma =  (1./(2.*pi_*infrared_Regulator*infrared_Regulator) -r/(2.*pi_*infrared_Regulator ) *gsl_sf_bessel_K1(r*infrared_Regulator))/(log(4./(infrared_Regulator*infrared_Regulator*r*r))) ; //Here we have to define the bessel function! ;
      
        
        number_type Qs2bar = Q2; //sqrt(coupling*coupling/casimir*OnePoint(b1,b2)*36.);
        number_type Qs2=pi_*8.*modifiedGamma*Qs2bar;
        //std::cout << "scale: " << scale << "\n";
        //std::cout << "Qs2: " << Qs2 << "\n";


		return (2./(pow(coupling,4.)*pow(r,8.))*(16.*exp(-Qs2)+32.*exp(-Qs2/2.)-64*exp(-3.*Qs2/4)-4.*exp(-1./4.*Qs2)*(Qs2*Qs2-2.*Qs2bar*Qs2bar*pow(r,4.)+8.*Qs2+48)+1./8.*exp(-Qs2/2)*(Qs2*Qs2*Qs2*Qs2+(4.*Qs2*Qs2+128.)*(2.*Qs2)+16.*(2.*Qs2)*(2.*Qs2)+1024.)+2.*(Qs2bar*Qs2bar*pow(r,4.)*(Qs2-4.)+40.)))/scale/scale;
        

	}

	


private:
	const gsl_interp2d_type* OnePoint_interpolator_;
	gsl_spline2d* OnePoint_spline_;
	gsl_interp_accel* xacc_;
	gsl_interp_accel* yacc_;
	number_type gridmax_;
	number_type gridstep_;

	number_type pi_;
	//
	number_type m_;
};


#endif
