#include <iostream>
#include "convolution.hpp"


int main ()
{
	typedef std::size_t size_type;
	typedef double number_type;

	Convolution<number_type> convolution(10.,0.2);

	convolution.get_convolved("profiles_averaged_20-21.txt");

	convolution.print_OnePoint("outfile.txt", 100);

	return 0;
}