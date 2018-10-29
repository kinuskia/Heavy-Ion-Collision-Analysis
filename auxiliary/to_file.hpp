#ifndef TO_FILE_HPP
#define TO_FILE_HPP

// Function which takes an STL vector of Vector-columns and writes them to a text file
template<typename number_type>
void to_file(std::string filename, std::vector<Vector<number_type>> data_columns)
{
	// using the Storage object from my numerical_analysis library
	Storage<number_type> data(data_columns.size());
	for (size_type i = 0; i < data_columns[0].size(); ++i)
	{
		for (size_type j = 0; j < data_columns.size(); ++j)
		{
			if (isnan((data_columns[j])[i]))
			{
				data.read_in(0);
			}
			else
			{
				data.read_in((data_columns[j])[i]);
			}
		}
	}
	data.write(filename, false);
}
// overloaded function for the special case of three columns
template<typename number_type>
void to_file(std::string filename, Vector<number_type> r, Vector<number_type> E, Vector<number_type> dE)
{
	std::vector<Vector<number_type>> columns(3);
	columns[0] = r;
	columns[1] = E;
	columns[2] = dE;
	to_file(filename, columns);
}

#endif