#include <iostream>
#include <vector>

#include "thickness.hpp"
#include "../auxiliary/to_file.hpp"


int main ()
{
	typedef std::size_t size_type;
	typedef double number_type;

	Thickness<number_type> thickness(0.25, 1.24, 9.604);

	// Read in W(r) function
 	thickness.initialize_W("weight_functions.txt");

 	// Normalize impact parameter distribution
	thickness.normalize();

	// save distribution in output file
	thickness.print_dist("b_dist.txt", 100);

	// Find percentiles of distribution
	std::vector<number_type> percentiles(101);
	thickness.find_percentiles(percentiles);

	// Save percentiles and the corresponding y-values in file
	std::vector<number_type> y_values(101);
	for (size_type i = 0; i < percentiles.size(); ++i)
	{
		y_values[i] = thickness.b_dist(percentiles[i]);
	}

	std::vector<std::vector<number_type>> output(2);
	output[0] = percentiles;
	output[1] = y_values;

	to_file("percentiles.txt", output);



	return 0;
}