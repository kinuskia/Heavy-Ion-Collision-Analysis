#ifndef TO_FILE_HPP
#define TO_FILE_HPP
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <vector>
#include <fstream>
#include <iomanip>
#include <cmath>

// Function which takes an GSL vector of Vector-columns and writes them to a text file

void to_file(std::string filename, std::vector<gsl_vector*> data_columns)
{
	typedef double number_type;
	typedef std::size_t size_type;
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

// overloaded function for the special case of two columns
void to_file(std::string filename, gsl_vector* r, gsl_vector* E)
{
	typedef double number_type;
	std::vector<gsl_vector*> columns(0);
	columns.push_back(r);
	columns.push_back(E);
	to_file(filename, columns);
}

// overload for STL vectors
template<typename number_type>
void to_file(std::string filename, std::vector<std::vector<number_type>> data_columns)
{
	typedef std::size_t size_type;
	std::ofstream outfile(filename);
	for (size_type i = 0; i < data_columns[0].size(); ++i)
	{
		for (size_type j = 0; j < data_columns.size(); ++j)
		{
			number_type value = (data_columns[j])[i] ;
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

// overload for gsl matrices
void to_file(std::string filename, gsl_matrix* data)
{
	typedef std::size_t size_type;
	typedef double number_type;
	std::vector<gsl_vector*> data_columns(0);
	for (size_type j = 0; j < data->size2; ++j)
	{
		gsl_vector* column = gsl_vector_alloc(data->size1);
		for (size_type i = 0; i < data->size1; ++i)
		{
			gsl_vector_set(column, i, gsl_matrix_get(data, i, j));
		}
		data_columns.push_back(column);
	}
	to_file(filename, data_columns);
	for (size_type i = 0; i < data_columns.size(); ++i)
	{
		gsl_vector_free(data_columns[i]);
	}
}

#endif