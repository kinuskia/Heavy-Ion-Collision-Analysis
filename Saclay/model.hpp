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
	// , gridmax_(10)
	// , gridstep_(0.2)
	{}

	// // initialize One-Point function
	// void initialize_OnePoint(std::string filename, size_type N, const gsl_interp2d_type* interp_method, number_type gridmax, number_type gridstep)
	// {
	// 	gsl_matrix* profile = gsl_matrix_alloc(N, N);
	// 	number_type dummy;
	// 	read_data(filename, profile, dummy, dummy, dummy);
		

	// 	const gsl_interp2d_type* interpolator = interp_method;
	// 	OnePoint_interpolator_ = interpolator;
	// 	number_type* za = new number_type[N*N];
	// 	gsl_spline2d* spline = gsl_spline2d_alloc(interpolator, N, N);
	// 	OnePoint_spline_ = spline;
	// 	gsl_interp_accel* xacc = gsl_interp_accel_alloc();
	// 	gsl_interp_accel* yacc = gsl_interp_accel_alloc();
	// 	xacc_ = xacc;
	// 	yacc_ = yacc;

	// 	// set interpolation grid values
	// 	number_type* x_sites = new number_type[N];
	// 	number_type* y_sites = new number_type[N];
	// 	gridmax_ = gridmax;
	// 	gridstep_ = gridstep;
	// 	number_type min = -gridmax_ +gridstep_*0.51 ;
	// 	number_type max = gridmax_ -gridstep_*0.51 ;
	// 	for (size_type i = 0; i < N ; ++i)
	// 	{
	// 		x_sites[i] = min + (max-min)*i/(N-1);
	// 		y_sites[i] = x_sites[i];
	// 	}

	// 	for (size_type i = 0; i < N; ++i)
	// 	{
	// 		for (size_type j = 0; j < N; ++j)
	// 		{
	// 			gsl_spline2d_set(OnePoint_spline_, za, i, j, gsl_matrix_get(profile, N-1 - j, i));
	// 		}
	// 	}

	// 	gsl_spline2d_init(OnePoint_spline_, x_sites, y_sites, za, N, N);

	// 	// free profile
	// 	gsl_matrix_free(profile);
	// }

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
		if (r > 9.604)
		{
			return 0;
		}
		return gsl_spline_eval(W_spline_, r, W_acc_);
	}

	number_type OnePoint(number_type R)
	{
		return W(R)/pi_;
	}

	// Define one-point position space correlation function
	// units in fm
	number_type OnePoint(number_type x, number_type y)
	{
		number_type R = sqrt(x*x+y*y);
		
		return OnePoint(R);
		// number_type result = gsl_spline2d_eval(OnePoint_spline_, x, y, xacc_, yacc_);
		// if (result < 0)
		// {
		// 	result = -result;
		// }

		// return result;

	}

	number_type TwoPoint(number_type r1, number_type r2, number_type b, number_type phiA, number_type phiB)
	{
		// // linear extrapolation if r near zero
		// number_type r_limit =0.02;
		// number_type result;
		
		// if (r < r_limit)
		// {
		// 	number_type r1 =0.021;
		// 	number_type r2 =0.03;
		// 	number_type y1 = TwoPoint_raw(R, r1);
		// 	number_type y2 = TwoPoint_raw(R, r2);

		// 	result = (y2-y1)/(r2-r1)*(r-r1)+y1;
		// }
		// else
		// {
		// 	result = TwoPoint_raw(R, r);
		// }

		// quadratic extrapolation if r near zero
		number_type r = sqrt(r1*r1 + r2*r2 - 2.*r1*r2*cos(phiA));
		number_type r_limit =0.02;
		number_type result;
		
		if (r < r_limit)
		{
			// number_type rA = 0.02;
			// number_type rB = 0.03;
			// number_type rC = 0.04;
			// number_type yA = TwoPoint_raw(R, rA);
			// number_type yB = TwoPoint_raw(R, rB);
			// number_type yC = TwoPoint_raw(R, rC);
			// number_type a = -(rC*yA-rC*yB-rB*yA+rB*yC-rA*yC+rA*yB)/(-rC+rB)/(-rC+rA)/(-rB+rA);
			// number_type b = -(-rA*rA*yC+rA*rA*yB+2.*rA*rB*yC-2.*rA*rB*yA+2.*rA*rC*yA-2.*rA*rC*yB+yB*rC*rC-rB*rB*yC-yA*rC*rC+rB*rB*yA)/(-rC+rB)/(-rC+rA)/(-rB+rA);
			// number_type c = yA;

			// result = a*(r-rA)*(r-rA) + b*(r-rA) + c;
			result = TwoPoint_raw(r1, r2, b, phiA, phiB);
		}
		else
		{
			result = TwoPoint_raw(r1, r2, b, phiA, phiB);
		}
		
		return result;	

	}

	// Two-point function with no extrapolation procedure for r->0
	/* naming of the angles
		phiA : phi(r1) - phi(r2)
		phiB : phi(r1) - phiR
	*/
	number_type TwoPoint_raw(number_type r1, number_type r2, number_type b, number_type phiA, number_type phiB)
	{
		number_type nc =3.;
        number_type hc = .1973;

        number_type r = sqrt(r1*r1 + r2*r2 - 2.*r1*r2*cos(phiA));

        number_type R = sqrt(r1*r1 + r2*r2 + 2.*r1*r2*cos(phiA));
        number_type RM = sqrt(R*R + b*b - 2.*b*(r1*cos(phiB) + r2*cos(phiA-phiB)));
        number_type RP = sqrt(R*R + b*b + 2.*b*(r1*cos(phiB) + r2*cos(phiA-phiB)));

     

        if (r > 6.) 
        {
        	return 0;
        } 
        r /= hc;

        // number_type R_cutoff = 0.01/hc;
        // if (r < R_cutoff)
        // {
        // 	r = R_cutoff;	
        // }

        // number_type b1= (x1+x2)/2.;
        // number_type b2= (y1+y2)/2.;
        number_type casimir= (nc*nc-1.)/2./nc;

        number_type alpha_s = 0.4095*0+0.25;
		number_type g = sqrt(4.*pi_*alpha_s);
		number_type W0 = W(0.);
		number_type Q0 = 1.1*0+1.24; //GeV
		number_type scale = Q0*Q0*Q0*Q0/3./alpha_s/W0;

		//number_type Q2 = Q0*Q0*sqrt(pi_*OnePoint(R/2)/W0);

		number_type Q2P = Q0*Q0*sqrt(pi_*OnePoint(RP/2)/W0);
		number_type Q2M = Q0*Q0*sqrt(pi_*OnePoint(RM/2)/W0);


       
        number_type coupling=g;
        number_type infrared_Regulator= m_; //0-1 1e-2
        //number_type modifiedGamma =  (1./(2.*pi_*infrared_Regulator*infrared_Regulator) -r/(2.*pi_*infrared_Regulator ) *gsl_sf_bessel_K1(r*infrared_Regulator))/(log(4./(infrared_Regulator*infrared_Regulator*r*r))) ; //Here we have to define the bessel function! ;
      
        
        //number_type Qs2bar = Q2; //sqrt(coupling*coupling/casimir*OnePoint(b1,b2)*36.);
       	number_type Qs2barP = Q2P;
       	number_type Qs2barM = Q2M;
        
        //number_type Qs2=pi_*8.*modifiedGamma*Qs2bar;
       	//number_type Qs2P = pi_*8.*modifiedGamma*Qs2barP;
       	//number_type Qs2M = pi_*8.*modifiedGamma*Qs2barM;

        number_type result0;
        // number_type result2;

        //result0  = 1./scale/scale*(1./(pow(coupling, 4.)*pow(r, 8.))*(16.*exp(-Qs2M/2-Qs2P/2) + 32.*exp(-Qs2P/2) - 64.*exp(-Qs2M/4-Qs2P/2) - 4.*exp(-Qs2P/4)*(Qs2P*Qs2P - 2.*Qs2barM*Qs2barM*pow(r, 4.) + 8.*Qs2P + 48.) + 1./8.*exp(-Qs2M/4-Qs2P/4)*(Qs2M*Qs2M*Qs2P*Qs2P + (4.*Qs2M*Qs2P + 128.)*(Qs2M + Qs2P) + 16.*(Qs2M+Qs2P)*(Qs2M+Qs2P) + 1024 )) + 2.*(Qs2barM*Qs2barM*pow(r, 4.)*(Qs2P - 4.) + 40.) );
        //result0 += 1./scale/scale*(1./(pow(coupling, 4.)*pow(r, 8.))*(16.*exp(-Qs2P/2-Qs2M/2) + 32.*exp(-Qs2M/2) - 64.*exp(-Qs2P/4-Qs2M/2) - 4.*exp(-Qs2M/4)*(Qs2M*Qs2M - 2.*Qs2barP*Qs2barP*pow(r, 4.) + 8.*Qs2M + 48.) + 1./8.*exp(-Qs2P/4-Qs2M/4)*(Qs2P*Qs2P*Qs2M*Qs2M + (4.*Qs2P*Qs2M + 128.)*(Qs2P + Qs2M) + 16.*(Qs2P+Qs2M)*(Qs2P+Qs2M) + 1024 )) + 2.*(Qs2barP*Qs2barP*pow(r, 4.)*(Qs2M - 4.) + 40.) );

		result0 = Large_Nc(r, infrared_Regulator, Qs2barM, Qs2barP, scale, coupling);


		//result2 = (2./(nc*nc*pow(coupling,4.)*pow(r,8.))*(2.*exp(-Qs2)*(2.*Qs2+8.)*(2.*Qs2+8.) + 4.*exp(-Qs2/2)*Qs2*(8.+Qs2) - 8.*exp(-3./4*Qs2)*(8.+Qs2)*(4.+Qs2) + 4.*exp(-Qs2/4)*(Qs2*Qs2-2.*Qs2bar*Qs2bar*pow(r,4.)+8.*Qs2+16.*Qs2) - 1./8*exp(-Qs2/2)*(Qs2*Qs2*Qs2*Qs2+(4.*Qs2*Qs2+128.)*(2.*Qs2)+16.*4.*Qs2*Qs2-1024.) - 2.*(Qs2bar*Qs2bar*pow(r,4.)*(Qs2-4.)+32.*Qs2-4.*Qs2*Qs2) ))/scale/scale;

		return result0;
	}

	// // Define connected (!) position space two-point correlation function
	// number_type TwoPoint(number_type x1, number_type y1, number_type x2, number_type y2)
	// {
 //     	number_type r = sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
 //     	number_type R = sqrt((x2+x1)*(x2+x1)+(y2+y1)*(y2+y1));

 //     	return TwoPoint(R, r);    

	// }

	number_type Large_Nc(number_type r, number_type infrared_Regulator, number_type Qs2barM, number_type Qs2barP, number_type scale, number_type coupling)
	{
		number_type hc = 0.1973;
		number_type r_limit = 0.02/hc;
		number_type result;
		if (r < r_limit)
		{
			number_type r1 = 0.02/hc;
			number_type r2 = 0.03/hc;
			number_type r3 = 0.04/hc;
			number_type y1 = Large_Nc_raw(r1, infrared_Regulator, Qs2barM, Qs2barP, scale, coupling);
			number_type y2 = Large_Nc_raw(r2, infrared_Regulator, Qs2barM, Qs2barP, scale, coupling);
			number_type y3 = Large_Nc_raw(r3, infrared_Regulator, Qs2barM, Qs2barP, scale, coupling);
			number_type a = -(r3*y1-r3*y2-r2*y1+r2*y3-r1*y3+r1*y2)/(-r3+r2)/(-r3+r1)/(-r2+r1);
			number_type b = -(-r1*r1*y3+r1*r1*y2+2.*r1*r2*y3-2.*r1*r2*y1+2.*r1*r3*y1-2.*r1*r3*y2+y2*r3*r3-r2*r2*y3-y1*r3*r3+r2*r2*y1)/(-r3+r2)/(-r3+r1)/(-r2+r1);
			number_type c = y1;

			result = a*(r-r1)*(r-r1) + b*(r-r1) + c;
		}
		else
		{
			result = Large_Nc_raw(r, infrared_Regulator, Qs2barM, Qs2barP, scale, coupling);
		}

		return result;
	}

	number_type Large_Nc_raw(number_type r, number_type infrared_Regulator, number_type Qs2barM, number_type Qs2barP, number_type scale, number_type coupling)
	{
		number_type modifiedGamma =  (1./(2.*pi_*infrared_Regulator*infrared_Regulator) -r/(2.*pi_*infrared_Regulator ) *gsl_sf_bessel_K1(r*infrared_Regulator))/(log(4./(infrared_Regulator*infrared_Regulator*r*r))) ;
		//number_type Qs2=pi_*8.*modifiedGamma*Qs2bar;
		number_type Qs2P = pi_*8.*modifiedGamma*Qs2barP;
       	number_type Qs2M = pi_*8.*modifiedGamma*Qs2barM;
		//number_type modifiedGamma =  (1./(2.*pi_*infrared_Regulator*infrared_Regulator) -r/(2.*pi_*infrared_Regulator ) *gsl_sf_bessel_K1(r*infrared_Regulator))/(log(4./(infrared_Regulator*infrared_Regulator*r*r))) ;
		number_type result0 = 0;
		number_type nc = 3;
		number_type Qs2 = Qs2P;
		number_type Qs2bar = Qs2barM;
		//result0  = (2./(pow(coupling,4.)*pow(r,8.))*(16.*exp(-Qs2)+32.*exp(-Qs2/2.)-64*exp(-3.*Qs2/4)-4.*exp(-1./4.*Qs2)*(Qs2*Qs2-2.*Qs2bar*Qs2bar*pow(r,4.)+8.*Qs2+48)+1./8.*exp(-Qs2/2)*(Qs2*Qs2*Qs2*Qs2+(4.*Qs2*Qs2+128.)*(2.*Qs2)+16.*(2.*Qs2)*(2.*Qs2)+1024.)+2.*(Qs2bar*Qs2bar*pow(r,4.)*(Qs2-4.)+40.)))/scale/scale;
		//result0 += (2./(nc*nc*pow(coupling,4.)*pow(r,8.))*(2.*exp(-Qs2)*(2.*Qs2+8.)*(2.*Qs2+8.) + 4.*exp(-Qs2/2)*Qs2*(8.+Qs2) - 8.*exp(-3./4*Qs2)*(8.+Qs2)*(4.+Qs2) + 4.*exp(-Qs2/4)*(Qs2*Qs2-2.*Qs2bar*Qs2bar*pow(r,4.)+8.*Qs2+16.*Qs2) - 1./8*exp(-Qs2/2)*(Qs2*Qs2*Qs2*Qs2+(4.*Qs2*Qs2+128.)*(2.*Qs2)+16.*4.*Qs2*Qs2-1024.) - 2.*(Qs2bar*Qs2bar*pow(r,4.)*(Qs2-4.)+32.*Qs2-4.*Qs2*Qs2) ))/scale/scale;
		
		result0  = 1./scale/scale*(1./(pow(coupling, 4.)*pow(r, 8.))*(16.*exp(-Qs2M/2-Qs2P/2) + 32.*exp(-Qs2P/2) - 64.*exp(-Qs2M/4-Qs2P/2) - 4.*exp(-Qs2P/4)*(Qs2P*Qs2P - 2.*Qs2barM*Qs2barM*pow(r, 4.) + 8.*Qs2P + 48.) + 1./8.*exp(-Qs2M/4-Qs2P/4)*(Qs2M*Qs2M*Qs2P*Qs2P + (4.*Qs2M*Qs2P + 128.)*(Qs2M + Qs2P) + 16.*(Qs2M+Qs2P)*(Qs2M+Qs2P) + 1024 ) + 2.*(Qs2barM*Qs2barM*pow(r, 4.)*(Qs2P - 4.) + 40.) ) );
        

        result0 += 1./scale/scale*(1./(pow(coupling, 4.)*pow(r, 8.))*(16.*exp(-Qs2P/2-Qs2M/2) + 32.*exp(-Qs2M/2) - 64.*exp(-Qs2P/4-Qs2M/2) - 4.*exp(-Qs2M/4)*(Qs2M*Qs2M - 2.*Qs2barP*Qs2barP*pow(r, 4.) + 8.*Qs2M + 48.) + 1./8.*exp(-Qs2P/4-Qs2M/4)*(Qs2P*Qs2P*Qs2M*Qs2M + (4.*Qs2P*Qs2M + 128.)*(Qs2P + Qs2M) + 16.*(Qs2P+Qs2M)*(Qs2P+Qs2M) + 1024 ) + 2.*(Qs2barP*Qs2barP*pow(r, 4.)*(Qs2M - 4.) + 40.) ) );


		return result0;
	}


private:
	const gsl_interp_type* W_interpolator_;
	gsl_spline* W_spline_;
	gsl_interp_accel* W_acc_;

	// const gsl_interp2d_type* OnePoint_interpolator_;
	// gsl_spline2d* OnePoint_spline_;
	// gsl_interp_accel* xacc_;
	// gsl_interp_accel* yacc_;
	// number_type gridmax_;
	// number_type gridstep_;

	number_type pi_;
	//
	number_type m_;
};


#endif
