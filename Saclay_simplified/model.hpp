#ifndef MODEL_HPP
#define MODEL_HPP

#include <iostream>
#include <cmath>
#include <gsl/gsl_spline.h>

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

		return W(r)/3.1415926;
	}

	// Define connected (!) position space two-point correlation function TwoPoint(x1, y1) = TwoPoint(x1, y1)*delta(x2-x1)*delta(y2-y1)
	number_type TwoPoint(number_type x, number_type y)
	{
		number_type r = sqrt(x*x+y*y);
		number_type Q2 = sqrt(OnePoint(x, y));

		return Q2*Q2*Q2*log(Q2/m_);
	}

private:
	const gsl_interp_type* W_interpolator_;
	gsl_spline* W_spline_;
	gsl_interp_accel* W_acc_;

	//
	number_type m_;
};


#endif