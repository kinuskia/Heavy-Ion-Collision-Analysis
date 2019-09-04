#ifndef READ_DATA_HPP
#define READ_DATA_HPP

#include <fstream>
#include <gsl/gsl_matrix.h>

typedef std::size_t size_type;
/* free functions to read in a data file and save it in a matrix */

template<typename number_type>
void insert(number_type number, size_type counter, gsl_matrix* data)
{
	size_type n_cols = data->size2;
	size_type rownumber = counter/n_cols;
	size_type colnumber = counter%n_cols;
	//assert(data->size1 == data->size2);
	
	gsl_matrix_set(data, rownumber, colnumber, number);
}

template<typename number_type>
void read_data(std::string filename, gsl_matrix* data, number_type & impact_parameter, number_type & n_participants, number_type & multiplicity)
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

		if (line[0] == '#') // handle comment lines
		{
			counter_lines++;
			if (line[2] == 'b') // get impact parameter
			{
				std::string impact_param_string = "";
				for (int i = 10; i < line.length(); ++i)
				{
					if (std::isalnum(line[i]) || (line[i] == '.') || (line[i] == '-'))
					{
						impact_param_string += line[i];
					}
				}
				impact_parameter = std::stod(impact_param_string);
				continue;
			}
			if (line[2] == 'n') // get number of wounded nucleons
			{
				std::string n_part_string = "";
				for (int i = 10; i < line.length(); ++i)
				{
					if (std::isalnum(line[i]) || (line[i] == '.') || (line[i] == '-'))
					{
						n_part_string += line[i];
					}
				}
				n_participants = std::stod(n_part_string);
				continue;
			}
			if (line[2] == 'm') // get multiplicity
			{
				std::string mult_string = "";
				for (int i = 10; i < line.length(); ++i)
				{
					if (std::isalnum(line[i]) || (line[i] == '.') || (line[i] == '-'))
					{
						mult_string += line[i];
					}
				}
				multiplicity = std::stod(mult_string);
				continue;
			}
			else
			{
				continue;
			}
		}
		

		bool end_of_number = true;
		std::string numberstring = "";
		bool last_char_is_number;
		for (int i = 0; i < line.length(); ++i)
		{
	
			// Detect beginning of a number
			bool is_number = std::isalnum(line[i]) || (line[i] == '.') || (line[i] == '-');
			//bool is_seperator = (line[i] == ' ' && line[i] == '\t');
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

// overload if one is just interested in the collision specs but not the profiles
template<typename number_type>
void read_data(std::string filename, number_type & impact_parameter, number_type & n_participants, number_type & multiplicity)
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

		if (line[0] == '#') // handle comment lines
		{
			counter_lines++;
			if (line[2] == 'b') // get impact parameter
			{
				std::string impact_param_string = "";
				for (int i = 10; i < line.length(); ++i)
				{
					if (std::isalnum(line[i]) || (line[i] == '.') || (line[i] == '-'))
					{
						impact_param_string += line[i];
					}
				}
				impact_parameter = std::stod(impact_param_string);
				continue;
			}
			if (line[2] == 'n') // get number of wounded nucleons
			{
				std::string n_part_string = "";
				for (int i = 10; i < line.length(); ++i)
				{
					if (std::isalnum(line[i]) || (line[i] == '.') || (line[i] == '-'))
					{
						n_part_string += line[i];
					}
				}
				n_participants = std::stod(n_part_string);
				continue;
			}
			if (line[2] == 'm') // get multiplicity
			{
				std::string mult_string = "";
				for (int i = 10; i < line.length(); ++i)
				{
					if (std::isalnum(line[i]) || (line[i] == '.') || (line[i] == '-'))
					{
						mult_string += line[i];
					}
				}
				multiplicity = std::stod(mult_string);
				continue;
			}
			else
			{
				continue;
			}
		}
		else
		{
			break;
		}
		
	}

}



#endif
