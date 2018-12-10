#ifndef BESSEL_DERIV_ZERO
#define BESSEL_DERIV_ZERO

#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_roots.h>
#include <assert.h>

// Free functions to compute the l'th positive zero crossing of Jm'(x)=0

struct bessel_deriv_params
{
	int m;
};

double bessel_deriv(double r, void * p)
{
	bessel_deriv_params* params = (bessel_deriv_params*) p;
	int m = params->m;

	return 1./2*(gsl_sf_bessel_Jn(m-1, r)-gsl_sf_bessel_Jn(m+1, r));
}


double find_bessel_deriv_root(int m, int l, double eps_rel, int max_iter)
{
	typedef double number_type;
	typedef std::size_t size_type;
	// solution and search region definitions
	number_type r;
	number_type lower;
	number_type upper;
	/* Using the zero crossings of Jm(x) as search region */
	if (m == 0)
	{
		lower = gsl_sf_bessel_zero_Jnu(0, l);
		upper = gsl_sf_bessel_zero_Jnu(0, l+1);
	}
	else if (l == 1)
	{
		lower = 0.1; // don't use lower = 0 because Jm'(0)=0 for m>1
		upper = gsl_sf_bessel_zero_Jnu(m, 1);
	}
	else 
	{
		lower = gsl_sf_bessel_zero_Jnu(m, l-1);
		upper = gsl_sf_bessel_zero_Jnu(m, l);
	}

	// root finding objects
	bessel_deriv_params params = {m};
	gsl_function F;
	F.function = &bessel_deriv;
	F.params = &params;

	const gsl_root_fsolver_type* T;
	gsl_root_fsolver* s;
	T = gsl_root_fsolver_brent;
	s = gsl_root_fsolver_alloc(T);
	gsl_root_fsolver_set(s, &F, lower, upper);


	// Do iterations

	for (size_type i = 0; i < max_iter; ++i)
	{
		int status = gsl_root_fsolver_iterate(s);
		r = gsl_root_fsolver_root(s);
		lower = gsl_root_fsolver_x_lower(s);
		upper = gsl_root_fsolver_x_upper(s);
		status = gsl_root_test_interval(lower, upper, 0, eps_rel); // estimate uncertainty
		number_type error = (upper-lower)/2.;
		//std::cout << "Iteration " << i << ": " << r << " +/- " << error << "\n"; 
		if (status == GSL_SUCCESS)
		{
			break;
		}
		assert(i != (max_iter-1)); // assert that iteration converges
	}

	gsl_root_fsolver_free(s);

	return r;

}


#endif