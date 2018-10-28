#include <iostream>
#include "/Users/Kianusch/Documents/Numerical_analysis/auxiliary_files/vector.hpp"
#include "/Users/Kianusch/Documents/Numerical_analysis/auxiliary_files/matrix.hpp"
#include "collision.hpp"
#include "polar.hpp"
#include <sstream>

#include <algorithm>

// function which converts string to size_t
template<typename stringtype>
std::size_t to_size_t(stringtype number)
{
	std::stringstream sstream(number);
	std::size_t result;
	sstream >> result;
	return result;
}

// get the k-th last character of a string
char last_character(std::string str, int k = 0)
{
	const int len = str.length();
	assert(len >=k);
	return str[len-1-k];
}



int main (int argc, char* argv[]) // command-line input: filename_begin, fileformat, impact_parameter, # of files
{
	typedef double number_type;
	typedef std::size_t size_type;

	Collision<number_type> PbPb(10, .2);

	std::string filename = argv[1];
	std::string impact_parameter = argv[3];
	filename += impact_parameter;
	filename += "/";
	std::string fileformat = argv[2];
	size_type n_files = to_size_t(argv[4]);
	// Read in Trento data
	PbPb.read_in(filename, fileformat, n_files);

	PbPb.normalize(1);
	PbPb.centralize();

	// compute phi-averaged energy density profile and save to text file
	Vector<number_type> radii(500);
	Vector<number_type> energy(radii.size());
	Vector<number_type> energy_err(radii.size());
	PbPb.average_azimuthal(radii, energy, energy_err);
	
	Storage<number_type> energies_averaged(3);
	for (size_type i = 0; i < radii.size(); ++i)
	{
		energies_averaged.read_in(radii[i]);
		energies_averaged.read_in(energy[i]);
		if (isnan(energy_err[i]))
		{
			energies_averaged.read_in(0);
		}
		else
		{
			energies_averaged.read_in(energy_err[i]);
		}
	}

	std::string outfile = "energies_averaged";
	outfile += impact_parameter;
	outfile += ".txt";
	energies_averaged.write(outfile,false);


	//PbPb.print(0);


	//std::cout << PbPb.interpolate(0, 0, 0) << "\n";


	return 0;
}