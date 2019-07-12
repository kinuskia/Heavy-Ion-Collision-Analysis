#ifndef MODEL_HPP
#define MODEL_HPP

#include <iostream>
#include <cmath>

#include <gsl/gsl_spline.h>
#include <gsl/gsl_sf_bessel.h>


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

	// getter for W(r) function
	number_type W(number_type r)
	{
		return gsl_spline_eval(W_spline_, r, W_acc_);
	}


	// Define one-point position space correlation function
	// units in fm
	number_type OnePoint(number_type x, number_type y)
	{
		number_type r = sqrt(x*x+y*y);
        //Here we should to the proper scale setting 
		return W(r)/3.1415926;
	}

	// Define connected (!) position space two-point correlation function
	number_type TwoPoint(number_type x1, number_type y1, number_type x2, number_type y2)
	{
        number_type nc =3.;
        number_type hc = .1973;
        number_type r = sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
        r /= hc;

        number_type R_cutoff = 0.1/hc;
        if (r < R_cutoff)
        {
        	r = R_cutoff;
        }

        number_type b1= (x1+x2)/2.;
        number_type b2= (y1+y2)/2.;
        number_type casimir= (nc*nc-1.)/2./nc;

        number_type alpha_s = 0.4095;
		number_type g = sqrt(4.*3.1415926*alpha_s);
		number_type W0 = W(0.);
		number_type Q0 = 1.1; //GeV
		number_type scale = Q0*Q0*Q0*Q0/3./alpha_s/W0;
		number_type Q2 = Q0*Q0*sqrt(3.1415926*OnePoint(b1,b2)/W0);
       
        number_type coupling=g;
        number_type infrared_Regulator= m_*1e-2;
        number_type modifiedGamma =  (1./(2.*3.1415926*infrared_Regulator*infrared_Regulator) -r/(2.*3.1415926*infrared_Regulator ) *gsl_sf_bessel_K1(r*infrared_Regulator))/(log(4./(infrared_Regulator*infrared_Regulator*r*r))) ; //Here we have to define the bessel function! ;
        number_type Qs2bar = Q2; //sqrt(coupling*coupling/casimir*OnePoint(b1,b2)*36.);
        number_type Qs2=3.1415926*8.*modifiedGamma*Qs2bar;
        //std::cout << "scale: " << scale << "\n";
        //std::cout << "Qs2: " << Qs2 << "\n";


		return (2./(pow(coupling,4.)*pow(r,8.))*(16.*exp(-Qs2)+32.*exp(-Qs2/2.)-64*exp(-3.*Qs2/4)-4.*exp(-1./4.*Qs2)*(Qs2-2.*Qs2bar*Qs2bar*pow(r,4.)+8.*Qs2+48)+1./8.*exp(-Qs2/2)*(Qs2*Qs2*Qs2*Qs2+(4.*Qs2*Qs2+128.)*(2.*Qs2)+16.*(2.*Qs2)*(2.*Qs2)+1024.)+2.*(Qs2bar*Qs2bar*pow(r,4.)*(Qs2-4.)+40.)))/scale/scale;
        




	}

private:
	const gsl_interp_type* W_interpolator_;
	gsl_spline* W_spline_;
	gsl_interp_accel* W_acc_;

	//
	number_type m_;
};


#endif
