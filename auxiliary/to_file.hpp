#ifndef TO_FILE_HPP
#define TO_FILE_HPP
#include <gsl/gsl_vector.h>
#include <vector>
#include <fstream>
#include <iomanip>

// Function which takes an GSL vector of Vector-columns and writes them to a text file

void to_file(std::string filename, std::vector<gsl_vector*> data_columns)
{
	typedef double number_type;
	std::ofstream outfile(filename);
	for (size_type i = 0; i < data_columns[0]->size; ++i)
	{
		for (size_type j = 0; j < data_columns.size(); ++j)
		{
			number_type value = gsl_vector_get(data_columns[j], i);
			if (isnan(value))
			{
				outfile << 0;
			}
			else 
			{
				outfile << std::setprecision(14) << value;
			}
			if (j < data_columns.size()-1)
			{
				outfile << " ";
			}
		}
		outfile << "\n";
	}	
}
// overloaded function for the special case of three columns
void to_file(std::string filename, gsl_vector* r, gsl_vector* E, gsl_vector* dE)
{
	typedef double number_type;
	std::vector<gsl_vector*> columns(0);
	columns.push_back(r);
	columns.push_back(E);
	columns.push_back(dE);
	to_file(filename, columns);
}

#endif