#include <iostream>
#include "/Users/Kianusch/Documents/Numerical_analysis/auxiliary_files/vector.hpp"
#include "/Users/Kianusch/Documents/Numerical_analysis/auxiliary_files/matrix.hpp"
#include "collision.hpp"
#include "polar.hpp"

#include <algorithm>



int main ()
{
	typedef double number_type;
	typedef std::size_t size_type;

	Collision<number_type> PbPb(10, 2);
	PbPb.read_in("../Trento-Example/PbPb/0.dat");

	PbPb.normalize(1);
	PbPb.centralize();
	PbPb.average_azimuthal();

	PbPb.print(0);


	std::cout << PbPb.interpolate(0, 0, 0) << "\n";


	return 0;
}