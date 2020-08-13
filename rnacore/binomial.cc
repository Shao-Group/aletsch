#include "binomial.h"

uint32_t compute_binomial_score(int n, double pr, int x)
{
	assert(x >= 0 && x <= n);
	if(x == 0) return UINT32_MAX;
	double p = compute_binomial_pvalue(n, pr, x);
	return (uint32_t)(-100.0 * log10(p));
}

double compute_binomial_pvalue(int n, double pr, int x)
{
	// compute the probability that observed >= x among n trials
	assert(x >= 0 && x <= n);
	if(x == 0) return 0;

	binomial_distribution<> b(n, pr);
	return cdf(complement(b, x - 1));
}
