#ifndef MODEL_HPP
#define MODEL_HPP

#include <iostream>
#include <cmath>

template<class REAL>
class Model
{
public:
	typedef std::size_t size_type;
	typedef REAL number_type;

	// constructor
	Model()
	{}

	// Define one-point position space correlation function
	// units in fm
	number_type OnePoint(number_type x, number_type y)
	{
		number_type r2 = x*x+y*y;
		number_type sigma = 2;

		// Gauss decay from origin
		return exp(-r2/(2*sigma*sigma));
	}

	// Define connected (!) position space two-point correlation function
	number_type TwoPoint(number_type x1, number_type y1, number_type x2, number_type y2)
	{
		number_type r = sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
		number_type s_c = 1.;
		size_type N = 170;

		number_type contact = exp(-r*r/(2*s_c*s_c))/2/3.1415926/s_c/s_c;
		return -OnePoint(x1, y1)*OnePoint(x2, y2)/N/(8*3.1415926)/(8*3.1415926) + 0.5*(OnePoint(x1, y1)+OnePoint(x2, y2))*contact/N/(8*3.1415926);
		//return OnePoint(x1, y1)*contact;

	}
};


#endif