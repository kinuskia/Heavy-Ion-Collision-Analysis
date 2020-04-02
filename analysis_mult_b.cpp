
#include "auxiliary/headers.hpp"

// Define a class to manage multiplicity distributions
template<typename REAL>
class MultDist
{
public:
	typedef REAL number_type;
	typedef std::size_t size_type;
	// constructor
	MultDist(std::vector<number_type> multiplicity, std::vector<number_type> probability)
	{

		// Construct spline
		n_bins_ = multiplicity.size();
		assert(probability.size() == n_bins_);
		mult_sites_ = new number_type[n_bins_+2];
		prob_sites_ = new number_type[n_bins_+2];

		for (size_type i = 0; i < n_bins_+2; ++i)
		{
			if (i == 0)
			{
				mult_sites_[0] = 3./2*multiplicity[0]-multiplicity[1]/2;
				prob_sites_[0] = 3./2*probability[0]-probability[1]/2;
			}
			else if (i == n_bins_+1)
			{
				mult_sites_[n_bins_+1] = 3./2*multiplicity[n_bins_-1]-multiplicity[n_bins_-2]/2;
				prob_sites_[n_bins_+1] = 3./2*probability[n_bins_-1]-probability[n_bins_-2]/2;
			}
			else
			{
				mult_sites_[i] = multiplicity[i-1];
				prob_sites_[i] = probability[i-1];
			}
		}


		const gsl_interp_type* interpolator = gsl_interp_cspline;
		interpolator_ = interpolator;
		gsl_spline* mult_prob_dist = gsl_spline_alloc(interpolator_, n_bins_+2);
		mult_prob_dist_ = mult_prob_dist;
		gsl_interp_accel* acc = gsl_interp_accel_alloc();
		acc_ = acc;
		gsl_spline_init(mult_prob_dist_, mult_sites_, prob_sites_, n_bins_+2);

		// initialize some helpful constants
		xmin_ = 3./2*multiplicity[0]-multiplicity[1]/2;
		xmax_ = 3./2*multiplicity[n_bins_-1]-multiplicity[n_bins_-2]/2;


	}

	number_type evaluate(number_type x)
	{
		return gsl_spline_eval(mult_prob_dist_, x, acc_);
	}

	// Integration step size is chosen such that new integral area dA remains constant
	number_type get_stepsize(number_type x, number_type dA)
	{
		number_type h = 0.01;
		// estimate (negative) derivate at x
		number_type D;
		if (x >= xmax_-h)
		{
			D = (evaluate(x)-evaluate(x-h))/h;
		}
		else if (x <= xmin_+h)
		{
			D = (evaluate(x+h)-evaluate(x))/h;
		}
		else
		{
			D = (evaluate(x+h)-evaluate(x-h))/2./h;
		}
		D = abs(D);

		number_type y = abs(evaluate(x));
		
		return 2.*dA*y/(2.*y*y-dA*D);
	}

	// integrate over probality distribution
	number_type integrate(number_type x_lower, number_type x_upper, size_type N)
	{
		assert(x_lower <= x_upper);
		assert(xmin_ <= x_lower);
		assert(x_upper <= xmax_);

		number_type integral = 0;
		for (size_type i = 0; i <= N; ++i)
		{
			number_type x = x_lower + (x_upper-x_lower)*i/N;
			number_type f = evaluate(x);
			if ((i == 0) || (i == N))
			{
				integral += f/2;
			}
			else
			{
				integral += f;
			}
		}
		integral *= (x_upper-x_lower)/N;
		
		return integral;
	}

	// integrate over mult times probality distribution -> mean/average
	number_type integrate_mean(number_type x_lower, number_type x_upper, size_type N)
	{
		assert(x_lower <= x_upper);
		assert(xmin_ <= x_lower);
		assert(x_upper <= xmax_);

		number_type integral = 0;
		for (size_type i = 0; i <= N; ++i)
		{
			number_type x = x_lower + (x_upper-x_lower)*i/N;
			number_type f = x*evaluate(x);
			if ((i == 0) || (i == N))
			{
				integral += f/2;
			}
			else
			{
				integral += f;
			}
		}
		integral *= (x_upper-x_lower)/N;
		
		return integral;
	}



	// Identify edges of centrality classes
	void get_edges(std::vector<number_type> & edges)
	{
		// size of edges determines centrality width
		size_type n_classes = edges.size()-1;
		number_type centrality_width = 100./n_classes;
		assert(n_classes >= 1);

		number_type current_integral = 0;
		number_type total = integrate(xmin_, xmax_, 500);
		size_type N = 50000;
		size_type counter_edges = 1;
		edges[n_classes] = xmax_;
		edges[0] = xmin_;

		number_type current_x = xmax_*.95;
		//std::cout << " integration bounds: " << xmin_ << " -> " << xmax_ << "\n";
		// fill rest of the edges by successively integrating from xmax to some decreasing x 
		for (size_type i = 0; i < N; ++i)
		{
			number_type dA = 0.0001;
			number_type h = get_stepsize(current_x, dA);
			//std::cout << h  << " " << evaluate(current_x)<< "\n";
			current_x -= h;
			
			if (current_x < xmin_)
			{
				break;
			}

			size_type n = 10 + (xmax_- current_x)/(xmax_- xmin_)*(500-10);
			number_type current_integral = integrate(current_x, xmax_*.99999, n)/total*100.;

			//std::cout << current_x << " " << current_integral * 100. << "\n";
			//std::cout << i  << " " << current_x << " " << current_integral << "\n";

			if (current_integral > centrality_width*counter_edges)
			{
				edges[n_classes-counter_edges] = current_x;
				counter_edges++;
			}

			

		}

		// for (size_type i = 0; i < edges.size(); ++i)
		// {
		// 	std::cout << 100-i << " " << edges[i] << "\n";
		// }
		
	}

	void get_normalizations(const std::vector<number_type> & edges,  std::vector<number_type> & normalizations)
	{
		assert(edges.size() == normalizations.size()+1);

		for (size_type i = 0; i < normalizations.size(); ++i)
		{
			size_type n = 10 + (edges[i+1]-edges[i])/(xmax_- xmin_)*(500-10);
			normalizations[i] = integrate_mean(edges[i], edges[i+1], n)/integrate(edges[i], edges[i+1], n);
		}

	}




private:
	size_type n_bins_;
	number_type* mult_sites_;
	number_type* prob_sites_;
	const gsl_interp_type* interpolator_;
	gsl_spline* mult_prob_dist_;
	gsl_interp_accel* acc_;

	number_type xmin_;
	number_type xmax_;
};




int main (int argc, char* argv[]) // command-line input:  filename_begin, fileformat, # of files, # histogram bins
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
	size_type n_bins = to_size_t(argv[4]);


	// Read in and pre-process Trento data
	Collision<number_type> PbPb(10, .2); // Create Collision object

	PbPb.read_in(filename, fileformat, n_files, 1., false); // read in Trento event files

	// Print file with the following columns: impact parameter, number of participants, multiplicity
	PbPb.collision_specs_to_file("output/collision_specs.txt");


	// Compute histogram data points for multiplicity and impact parameter

	// Multiplicity
	std::vector<number_type> multiplicity(n_bins);
	std::vector<number_type> probability(n_bins);
	std::vector<number_type> probability_error(n_bins);
	PbPb.histogram_mult("output/mult_hist.txt", n_bins, multiplicity, probability, probability_error, true); // normed = true/false
	
	// Initialize distribution manager
	MultDist<number_type> multdist(multiplicity, probability);

	// Identify centrality intervals
	std::vector<number_type> edges(101);
	multdist.get_edges(edges);

	// get mean normalization of each centrality class
	std::vector<number_type> normalizations(edges.size()-1);
	multdist.get_normalizations(edges, normalizations);

	for (size_type i = 0; i < normalizations.size(); ++i)
	{
		std::cout << i << " : " << edges[normalizations.size()-i-1] << " -> " << edges[normalizations.size()-i] << " : " << normalizations[normalizations.size()-1-i] << "\n";
	}
	std::reverse(normalizations.begin(), normalizations.end());

	std::vector<std::vector<number_type>> output_norms(1);
	output_norms[0] = normalizations;
	to_file("output/average_multiplicity.txt", output_norms);







	//impact parameter
	std::vector<number_type> impact_params(n_bins);
	std::vector<number_type> prob_impact(n_bins);
	std::vector<number_type> prob_error_impact(n_bins);
	PbPb.histogram_b("output/b_hist.txt", n_bins , impact_params, prob_impact, prob_error_impact, true); // normed = true/false
	// Initialize distribution manager
	MultDist<number_type> impactdist(impact_params, prob_impact);
	// Identify centrality intervals
	std::vector<number_type> edges_impact(101);
	std::vector<number_type> edges_prob(101);
	impactdist.get_edges(edges_impact);

	for (size_type i = 0; i < edges_prob.size(); ++i)
	{
		edges_prob[i] = impactdist.evaluate(edges_impact[i]);
	}

	std::vector<std::vector<number_type>> output_impact(2);
	output_impact[0] = edges_impact;
	output_impact[1] = edges_prob;

	to_file("output/percentiles_b.txt", output_impact);

	return 0;
}