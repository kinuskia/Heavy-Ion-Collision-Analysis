#ifndef STATISTICS_HPP
#define STATISTICS_HPP

#include <gsl/gsl_vector.h>
#include <assert.h>
#include <cmath>

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

#endif