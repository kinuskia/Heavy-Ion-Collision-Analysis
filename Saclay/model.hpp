#ifndef MODEL_HPP
#define MODEL_HPP

#include <iostream>
#include <cmath>

#include <gsl/gsl_matrix.h>

#include <gsl/gsl_spline.h>
#include <gsl/gsl_sf_bessel.h>
#include "../auxiliary/bessel_deriv_zero.hpp"
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
	{}

	// initialize W(r) function
	void initialize_W(std::string filename)
	{
		number_type* W_sites;
		number_type* r_sites;
		// Count number of sites in text file
		std::ifstream counter_sites(filename);
		size_type Nr = 0;
		while (counter_sites)
		{
			std::string line;
			std::getline(counter_sites, line);
			if (line == "")
			{
				continue;  // ignore empty lines
			}
			Nr++;
		}
		counter_sites.close();

		W_sites = new number_type[Nr];
		r_sites = new number_type[Nr];

		// Read in W data from text file
		std::ifstream data(filename);
		std::string line;
		size_type i = 0;
		while (data)
		{
			std::getline(data, line); // Read in current line
			if (line == "")
			{
				continue;  // ignore empty lines
			}
			std::string str_radius = "";
			std::string str_W = "";
			bool radius_read_in = false;
			for (int i = 0; i < line.length(); ++i)
			{
				if (line[i] == ' ')
				{
					radius_read_in = true;
					continue;
				}
				if (!radius_read_in)
				{
					str_radius += line[i];
				}
				else
				{
					str_W += line[i];
				}
			}
			r_sites[i] = std::stod(str_radius);
			W_sites[i] = std::stod(str_W);
			++i;
		}
		data.close();

		//  Generate a spline
	
		const gsl_interp_type* interpolator = gsl_interp_cspline;
		gsl_spline* spline = gsl_spline_alloc(interpolator, Nr);
		gsl_interp_accel* acc = gsl_interp_accel_alloc();

		W_interpolator_ = interpolator;
		W_spline_ = spline;
		W_acc_ = acc;
		gsl_spline_init(W_spline_, r_sites, W_sites, Nr);
	}

	// import background coefficients
	void initialize_coeffs(std::string filename, size_type mMax, size_type lMax)
	{
		number_type dummy;
		gsl_matrix* data = gsl_matrix_alloc(mMax+1, lMax);
		read_data(filename, data, dummy, dummy, dummy);
		e_ml_ = data;
	}

	// getter for background coefficients
	number_type e_ml(size_type m, size_type l)
	{
		return gsl_matrix_get(e_ml_, m, l-1);
	}


	// getter for W(r) function
	number_type W(number_type r)
	{
		return gsl_spline_eval(W_spline_, r, W_acc_);
	}

	// integrate over W
	number_type integrate_W_r(number_type rMin, number_type rMax, size_type N)
	{
		number_type result = 0;
		number_type width = (rMax-rMin)/(N-1);
		for (size_type i = 0; i <= N; ++i)
		{
			number_type r = rMin + (rMax-rMin)*i/N;
			number_type f = W(r)*r;
			if ( (i==0) || (i==N) )
			{
				result += f/2;
			}
			else
			{
				result += f;
			}
		}
		result *= width;
		return result;
	}

	// initialize rho(r) function
	void initialize_rho()
	{
		number_type* rho_sites;
		number_type* r_sites;
		size_type Nr = 100;
		number_type rMax = 9.604;

		rho_sites = new number_type[Nr];
		r_sites = new number_type[Nr];

		for (size_type i = 0; i < Nr; ++i)
		{
			r_sites[i] = rMax*i/(Nr-1);
			size_type N = 5 + 195*i/Nr;
			rho_sites[i] = sqrt(2.*integrate_W_r(0,r_sites[i], N));
		}

		//  Generate a spline
	
		const gsl_interp_type* interpolator = gsl_interp_cspline;
		gsl_spline* spline = gsl_spline_alloc(interpolator, Nr);
		gsl_interp_accel* acc = gsl_interp_accel_alloc();

		rho_interpolator_ = interpolator;
		rho_spline_ = spline;
		rho_acc_ = acc;
		gsl_spline_init(rho_spline_, r_sites, rho_sites, Nr);
	}

	// getter for rho
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


	// Define one-point position space correlation function
	// units in fm
	number_type OnePoint(number_type x, number_type y)
	{
		number_type r = sqrt(x*x+y*y);
		number_type phi = atan2(y,x);
        number_type result = W(r)/3.1415926;

        size_type mMax = 10;
        size_type lMax = 10;
        for (size_type m=2; m<=mMax; ++m)
        {
        	// if (m%2 == 1)
        	// {
        	// 	continue;
        	// }
        	for (size_type l=1; l <= lMax; ++l)
        	{
        		result += 2.*W(r)*cos(phi*m)*Bessel(m, l, r)*e_ml(m, l);
        	}
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
		number_type g = sqrt(4.*3.1415926*alpha_s);
		number_type W0 = W(0.);
		number_type Q0 = 1.1*0+1.24; //GeV
		number_type scale = Q0*Q0*Q0*Q0/3./alpha_s/W0;
		number_type Q2 = Q0*Q0*sqrt(3.1415926*OnePoint(b1,b2)/W0);
       
        number_type coupling=g;
        number_type infrared_Regulator= m_*1e-2; //0-1 1e-2
        number_type modifiedGamma =  (1./(2.*3.1415926*infrared_Regulator*infrared_Regulator) -r/(2.*3.1415926*infrared_Regulator ) *gsl_sf_bessel_K1(r*infrared_Regulator))/(log(4./(infrared_Regulator*infrared_Regulator*r*r))) ; //Here we have to define the bessel function! ;
      
        
        number_type Qs2bar = Q2; //sqrt(coupling*coupling/casimir*OnePoint(b1,b2)*36.);
        number_type Qs2=3.1415926*8.*modifiedGamma*Qs2bar;
        //std::cout << "scale: " << scale << "\n";
        //std::cout << "Qs2: " << Qs2 << "\n";


		return (2./(pow(coupling,4.)*pow(r,8.))*(16.*exp(-Qs2)+32.*exp(-Qs2/2.)-64*exp(-3.*Qs2/4)-4.*exp(-1./4.*Qs2)*(Qs2*Qs2-2.*Qs2bar*Qs2bar*pow(r,4.)+8.*Qs2+48)+1./8.*exp(-Qs2/2)*(Qs2*Qs2*Qs2*Qs2+(4.*Qs2*Qs2+128.)*(2.*Qs2)+16.*(2.*Qs2)*(2.*Qs2)+1024.)+2.*(Qs2bar*Qs2bar*pow(r,4.)*(Qs2-4.)+40.)))/scale/scale;
        

	}

	


private:
	const gsl_interp_type* W_interpolator_;
	gsl_spline* W_spline_;
	gsl_interp_accel* W_acc_;

	const gsl_interp_type* rho_interpolator_;
	gsl_spline* rho_spline_;
	gsl_interp_accel* rho_acc_;

	gsl_matrix* bessel_deriv_zeros_;

	gsl_matrix* e_ml_;

	//
	number_type m_;
};


#endif
