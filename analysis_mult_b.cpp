
#include "auxiliary/headers.hpp"



int main (int argc, char* argv[]) // command-line input: filename_begin, fileformat, # of files
{
	typedef double number_type;
	typedef std::size_t size_type;

	// set timer
	std::time_t start = std::time(nullptr);

	// Evaluate command-line input
	std::string filename = argv[1];
	filename += "/";
	std::string fileformat = argv[2];
	size_type n_files = to_size_t(argv[3]);


	// Read in and pre-process Trento data
	Collision<number_type> PbPb(10, .2); // Create Collision object

	PbPb.read_in(filename, fileformat, n_files, false); // read in Trento event files

	// Print file with the following columns: impact parameter, number of participants, multiplicity
	PbPb.collision_specs_to_file("output/collision_specs.txt");


	// Compute histogram data points for multiplicity and impact parameter
	size_type N_bins = 50;
	PbPb.histogram_mult("output/mult_hist.txt", N_bins, true); // normed = true/false
	PbPb.histogram_b("output/b_hist.txt", N_bins, true); // normed = true/false

	
	return 0;
}