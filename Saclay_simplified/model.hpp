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

	// Define connected (!) position space two-point correlation function TwoPoint(x1, y1) = TwoPoint(x1, y1)*delta(x2-x1)*delta(y2-y1)
	number_type TwoPoint(number_type x, number_type y)
	{
		size_type N = 170;
		return pow(OnePoint(x, y),3./3)/N/(8*3.1415926);
	}
};


#endif