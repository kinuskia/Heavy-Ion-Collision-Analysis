#ifndef THICKNESS_HPP
#define THICKNESS_HPP

#include <gsl/gsl_spline.h>
#include "../auxiliary/read_data.hpp"
#include "../auxiliary/to_file.hpp"
#include "../auxiliary/to_size_t.hpp"


template <class REAL>
class Thickness
{
public:
	typedef std::size_t size_type;
	typedef REAL number_type;

	// // constructor
	Thickness(number_type alpha_s, number_type Q0, number_type rMax)
	: pi_(3.1415926)
	, alpha_s_(alpha_s)
	, Q0_(Q0)
	, rMax_(rMax)
	, normalization_(1.)
	, Nr_(25)
	, Nphi_(20)
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
			if (line[0] == '#') // ignore header
			{
				continue;
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
		if (r > rMax_)
		{
			return 0;
		}

		return gsl_spline_eval(W_spline_, r, W_acc_);
	}

	// getter for saturation scale
	number_type Q2(number_type r)
	{
		number_type g = sqrt(4.*pi_*alpha_s_);
		number_type W0 = W(0.);
		number_type scale = Q0_*Q0_*Q0_*Q0_/3./alpha_s_/W0;
		return sqrt(3.*g*g/4.*scale*W(r)/pi_);
	}

	// getter for one_point function in polar coords, d_phi=phi-phiR
	number_type OnePoint_r_phi(number_type r, number_type d_phi, number_type b)
	{
		number_type g = sqrt(4.*pi_*alpha_s_);
		number_type r_minus = sqrt(r*r + b*b/4. - b*r*cos(d_phi));
		number_type r_plus = sqrt(r*r + b*b/4. + b*r*cos(d_phi));
		return 4./3./g/g*Q2(r_minus)*Q2(r_plus);
	}

	// getter for OnePoint function, phi-dependence is integrated out
	number_type OnePoint_r(number_type r, number_type b)
	{
		// integrate OnePoint_phi over d_Phi with trapezoidal rule
		number_type result = 0;
		size_type N = Nphi_;
		number_type width = pi_/N;

		for (size_type i = 0; i <= N; ++i)
		{
			number_type phi = pi_*i/N;
			number_type f = OnePoint_r_phi(r, phi, b);
			if ((i == 0) || (i == N))
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

	// getter for One point function with r integrated out
	number_type OnePoint(number_type b)
	{
		// trapezoidal rule
		number_type result = 0; 
		size_type N = Nr_;
		number_type width = rMax_/N;

		for (size_type i = 0; i <= N; ++i)
		{
			number_type r = rMax_*i/N;
			number_type f = r*pow(OnePoint_r(r, b), 3./4);
			if ((i == 0) || (i == N))
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

	// get normalization constant for impact param distribution
	number_type normalization()
	{
		// trapezoidal rule
		number_type result = 0; 
		size_type N = 2.*Nr_;
		number_type width = 2.*rMax_/N;

		for (size_type i = 0; i <= N; ++i)
		{
			number_type b = 2.*rMax_*i/N;
			number_type f = OnePoint(b);
			if ((i == 0) || (i == N))
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

	// Compute normalization constant once and for all
	void normalize()
	{
		normalization_ = normalization();
	}

	// getter for impact param distribution
	number_type b_dist(number_type b)
	{

		return OnePoint(b)/normalization_;
	}

	// integrate OnePoint function over certain range of impact parameters
	number_type integ_OnePoint(number_type bmin, number_type bmax)
	{
		// trapezoidal rule
		number_type result = 0; 
		size_type N = Nphi_;
		number_type width = (bmax-bmin)/N;

		for (size_type i = 0; i <= N; ++i)
		{
			number_type b = bmin + (bmax-bmin)*i/N;
			number_type f = b_dist(b);
			if ((i == 0) || (i == N))
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

	// Newton's method of finding percentile intervals
	number_type find_next_bound(number_type b0)
	{
		size_type N = 10;
		number_type result = b0;
		for (size_type i = 0; i < N; ++i)
		{
			result -= (integ_OnePoint(b0, result)-1./100)/b_dist(result);
		}

		return result;
	}

	// find percentile bounds for b
	void find_percentiles (std::vector<number_type> & percentiles)
	{
		assert(percentiles.size()==101);
		percentiles[0] = 0;
		percentiles[100] = 2.*rMax_;
		for (size_type i = 1; i < 100; ++i)
		{
			percentiles[i] = find_next_bound(percentiles[i-1]);
		}

	}


	// print impact parameter dist to output file
	void print_dist(std::string filename, size_type N)
	{
		std::vector<number_type> b(N);
		std::vector<number_type> p(N);

		for (size_type i = 0; i < b.size(); ++i)
		{
			number_type b_current = 2.*rMax_*i/N;
			b[i] = b_current;
			p[i] = b_dist(b_current);
		}

		std::vector<std::vector<number_type>> output(2);
		output[0] = b;
		output[1] = p;

		to_file(filename, output);
	}

private:
	const gsl_interp_type* W_interpolator_;
	gsl_spline* W_spline_;
	gsl_interp_accel* W_acc_;

	number_type pi_;
	number_type alpha_s_;
	number_type Q0_;

	number_type rMax_;
	number_type normalization_;
	number_type Nr_;
	number_type Nphi_;





};

#endif