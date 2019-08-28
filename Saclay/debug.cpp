#include <iostream>

#include "model.hpp"

#include <vector>
#include <string>
#include "../auxiliary/to_file.hpp"

int main ()
{
	typedef double number_type;
	typedef std::size_t size_type;
	

	size_type N = 500;
	number_type r_min = 0.0001;
	number_type r_max = 2;
	size_type centrality_min = size_type(0);
	size_type centrality_max = size_type(1);
	number_type phi = 10*3.1415926/180;

	

	std::vector<number_type> r(N);
	std::vector<number_type> y(N);

	Model<number_type> LargeNc(0.14);
	std::string centrality = std::to_string(centrality_min) + "-" + std::to_string(centrality_max);

	LargeNc.initialize_W("weight_functions_"+centrality+".txt");
	
	for (size_type i = 0; i < N; ++i)
	{
		r[i] = r_min + (r_max-r_min)*i/N;
		y[i] = LargeNc.TwoPoint_test(r[i]*cos(phi), r[i]*sin(phi), 0, 0);
	}

	std::vector<std::vector<number_type>> output(2);
	output[0] = r;
	output[1] = y;

	to_file("TwoPoint_0-1.txt", output);




	return 0;
}