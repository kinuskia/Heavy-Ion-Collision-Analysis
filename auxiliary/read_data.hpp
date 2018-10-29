#ifndef READ_DATA_HPP
#define READ_DATA_HPP

#include <fstream>
#include "/Users/Kianusch/Documents/Numerical_analysis/auxiliary_files/vector.hpp"
#include "/Users/Kianusch/Documents/Numerical_analysis/auxiliary_files/matrix.hpp"

typedef std::size_t size_type;
/* free functions to read in a data file and save it in a matrix */

template<typename number_type>
void insert(number_type number, size_type counter, Matrix<number_type> & data)
{
	size_type n_cols = data.colsize();
	size_type rownumber = counter/n_cols;
	size_type colnumber = counter%n_cols;
	assert(data.colsize() == data.rowsize());
	
	data(rownumber, colnumber) = number;
}

template<typename number_type>
void read_data(std::string filename, Matrix<number_type> & data, size_type skip_rows = 0)
{
	std::ifstream infile(filename);
	std::string line;

	size_type counter_numbers = 0;
	size_type counter_lines = 0;
	while (infile)
	{
		std::getline(infile, line); // Read in current line
		if (line == "")
		{
			continue;  // ignore empty lines
		}

		if (counter_lines < skip_rows) // skip first rows if one wished so
		{
			counter_lines++;
			continue;
		}

		bool end_of_number = true;
		std::string numberstring = "";
		bool last_char_is_number;
		for (int i = 0; i < line.length(); ++i)
		{
			// Detect beginning of a number
			bool is_number = std::isalnum(line[i]) || (line[i] == '.') || (line[i] == '-');
			bool is_seperator = (line[i] == ' ' && line[i] == '\t');
			if (is_number) // Detect beginning of a number
			{
				end_of_number = false;
				numberstring += line[i];
			}
			else if (!end_of_number)
			{
				end_of_number = true;
				number_type number = std::stod(numberstring);
				insert(number, counter_numbers, data); // save number in the correct vector
				counter_numbers++;
				numberstring = "";
			}
			else
			{
				continue; // ignore second, third ... space in a series of spaces
			}
			if (i == line.length()-1)
			{
				last_char_is_number = is_number;
			}
		}
		if (last_char_is_number) // if data line does not end with a weird character
		{
			number_type number = std::stod(numberstring); // end of line finishes number
			insert(number, counter_numbers, data);
			counter_numbers++;
			numberstring = "";
		}
		counter_lines++;
	}

}

#endif
