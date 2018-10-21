#include <iostream>
#include "/Users/Kianusch/Documents/Numerical_analysis/auxiliary_files/vector.hpp"
#include "/Users/Kianusch/Documents/Numerical_analysis/auxiliary_files/matrix.hpp"
#include "collision.hpp"

int main ()
{
	typedef double number_type;
	typedef std::size_t size_type;

	Collision<number_type> PbPb(10, 0.2);
	PbPb.read_in("../Trento-Example/PbPb/0.dat");

	PbPb.normalize(1);
	PbPb.centralize();
	PbPb.centralize();

	//PbPb.print(0);


	return 0;
}