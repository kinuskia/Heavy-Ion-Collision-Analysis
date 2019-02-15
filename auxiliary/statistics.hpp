#ifndef STATISTICS_HPP
#define STATISTICS_HPP

#include <gsl/gsl_vector.h>
#include <assert.h>
#include <cmath>
#include <vector>
#include "complex_matrix.hpp"

// Compute mean
template<typename number_type>
void mean(gsl_vector* v, number_type & result)
{
	typedef std::size_t size_type;
	number_type sum = 0;
	size_type N = v->size;
	assert(N>0);
	for (size_type i = 0; i < N; ++i)
	{
		sum += gsl_vector_get(v, i);
	}

	result = sum / N;
}

// Compute mean (overload for STL vector)
template<typename number_type>
void mean(std::vector<number_type> v, number_type & result)
{
	typedef std::size_t size_type;
	number_type sum = 0;
	size_type N = v.size();
	assert(N>0);
	for (size_type i = 0; i < N; ++i)
	{
		sum += v[i];
	}

	result = sum / N;
}

// compute the (unbiased) standard deviation of a gsl vector
template<typename number_type>
void std_deviation(gsl_vector* v, number_type & result)
{
	typedef std::size_t size_type;
	size_type N = v->size;
	
	// Compute mean
	number_type v_mean;
	mean(v, v_mean);

	// compute standard deviation
	number_type sum_squared = 0;
	
	assert(N>1);
	for (size_type i = 0; i < N; ++i)
	{
		sum_squared += (gsl_vector_get(v, i)-v_mean)*(gsl_vector_get(v, i)-v_mean);
	}

	result = sqrt(sum_squared/(N-1));
}

// compute the (unbiased) standard deviation of an STL vector
template<typename number_type>
void std_deviation(std::vector<number_type> v, number_type & result)
{
	typedef std::size_t size_type;
	size_type N = v.size();
	
	// Compute mean
	number_type v_mean;
	mean(v, v_mean);

	// compute standard deviation
	number_type sum_squared = 0;
	
	assert(N>1);
	for (size_type i = 0; i < N; ++i)
	{
		sum_squared += (v[i]-v_mean)*(v[i]-v_mean);
	}

	result = sqrt(sum_squared/(N-1));
}

// compute the mean of a gsl vector as well as its uncertainty
template<typename number_type>
void mean(gsl_vector* v, number_type & result, number_type & result_error)
{
	typedef std::size_t size_type;
	size_type N = v->size;
	
	// compute mean
	mean(v, result);

	//compute uncertainty
	number_type v_std;
	std_deviation(v, v_std);

	result_error = v_std/sqrt(N);
}

// compute the mean of an STL vector as well as its uncertainty
template<typename number_type>
void mean(std::vector<number_type> v, number_type & result, number_type & result_error)
{
	typedef std::size_t size_type;
	size_type N = v.size();
	
	// compute mean
	mean(v, result);

	//compute uncertainty
	number_type v_std;
	std_deviation(v, v_std);

	result_error = v_std/sqrt(N);
}



// compute the mean of a STL vector of complex matrices as well as its uncertainty
template<typename number_type>
void mean (const std::vector<complex_matrix<number_type>> & data, complex_matrix<number_type> & matrix_mean, complex_matrix<number_type> & matrix_err)
{
	typedef std::size_t size_type;
	assert(data[0].colsize() == matrix_mean.colsize());
	assert(data[0].rowsize() == matrix_mean.rowsize());
	assert(data[0].colsize() == matrix_err.colsize());
	assert(data[0].rowsize() == matrix_err.rowsize());

	// compute mean
	matrix_mean = 0;
	for (size_type k = 0; k < data.size(); ++k)
	{
		matrix_mean += data[k];
	}
	matrix_mean /= data.size();

	// compute variance
	matrix_err = 0;
	for (size_type i = 0; i < matrix_err.rowsize(); ++i)
	{
		for (size_type j = 0; j < matrix_err.colsize(); ++j)
		{
			for (size_type k = 0; k < data.size(); ++k)
			{
			// variance of real part
				number_type diff_real = data[k].get_real(i, j) - matrix_mean.get_real(i, j);
			// variance of imaginary part
				number_type diff_imag = data[k].get_imag(i, j) - matrix_mean.get_imag(i, j);

				matrix_err.increase_entry(i, j, diff_real*diff_real + diff_imag*diff_imag, 0);
			}
		}
	}
	matrix_err /= (data.size() - 1);

	// Transform matrix of variances in errors
	matrix_err /= data.size();
	for (size_type i = 0; i < matrix_err.rowsize(); ++i)
	{
		for (size_type j = 0; j < matrix_err.colsize(); ++j)
		{
			number_type square = matrix_err.get_real(i, j);
			matrix_err.set_entry(i, j, sqrt(square), 0);
		}
	}


}

// signum function: returns +1 or -1 (input type is used as output)
template<typename number_type>
number_type sign(number_type value)
{
	number_type sign_value;
	if (value >= 0)
	{
		sign_value = number_type(1.0);
	}
	else
	{
		sign_value = number_type(-1.0);
	}
	return sign_value;
}

#endif