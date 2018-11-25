#include <iostream>
#include <fstream>
#include <string>

void generate_trento_script(std::string impact_parameter, std::string filename)
{
	// Create correct filename
		std::ofstream outfile(filename);

		// Specify the projectile option twice
		outfile << "projectile = Pb" << "\n";
		outfile << "projectile = Pb" << "\n";
		outfile << "number-events = 1000" << "\n";

		// don't print event properties to stdout, save to text file
		outfile << "quiet = true" << "\n";
		outfile << "output = /Users/Kianusch/Documents/Studium/Semester/WiSe1819/Bachelor-Arbeit/Heavy-Ion-Collision-Analysis/Trento/PbPb";
		outfile << impact_parameter << "\n";

		outfile << "reduced-thickness = 0" << "\n";
		outfile << "fluctuation = 1" << "\n";
		outfile << "nucleon-width = 0.5" << "\n";
		outfile << "nucleon-min-dist = 0" << "\n";
		outfile << "cross-section = 6.4" << "\n";
		outfile << "normalization = 1" << "\n";

		// impact parameter, leave commented out for min-bias
		outfile << "b-min = " << impact_parameter << "\n";
		outfile << "b-max = " << impact_parameter << "\n";

		outfile << "grid-max = 10" << "\n";
		outfile << "grid-step = 0.2" << "\n";

}

// Create a Trento collision script of given name and impact parameter

int main (int argc, char* argv[]) // command-line input: impact parameter, file name
{
	generate_trento_script(argv[1], argv[2]);

	return 0;
}