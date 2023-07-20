/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#include "constants.h"
#include "region.h"
#include "config.h"
#include "binomial.h"
#include <algorithm>

using namespace std;

/*
region::region(int32_t _lpos, int32_t _rpos, int _ltype, int _rtype)
	:lpos(_lpos), rpos(_rpos), mmap(NULL), imap(NULL), ltype(_ltype), rtype(_rtype), subregion_gap(-1), subregion_length(-1), subregion_overlap(-1)
{
} 
*/

region::region(int32_t _lpos, int32_t _rpos, int _ltype, int _rtype, const split_interval_map *_mmap, const split_interval_map *_imap, const parameters &c, const sample_profile &s)
	:lpos(_lpos), rpos(_rpos), mmap(_mmap), imap(_imap), ltype(_ltype), rtype(_rtype), cfg(c), sp(s)
{
	build_join_interval_map();
	if(ltype == RIGHT_SPLICE && rtype == LEFT_SPLICE) smooth_join_interval_map();
	split_large_region();
	build_partial_exons();
	//calculate_significance();
} 

region::~region()
{}

int region::build_join_interval_map()
{
	jmap.clear();

	PSIMI pei = locate_boundary_iterators(*mmap, lpos, rpos);
	SIMI lit = pei.first, rit = pei.second;

	if(lit == mmap->end() || rit == mmap->end()) return 0;

	SIMI it = lit;
	while(true)
	{
		//if(it->second >= 2) 
		jmap += make_pair(it->first, 1);
		if(it == rit) break;
		it++;
	}

	for(JIMI it = jmap.begin(); it != jmap.end(); it++)
	{
		assert(it->second == 1);
	}

	return 0;
}

int region::split_large_region()
{
	if(boost::distance(jmap.begin(), jmap.end()) != 1) return 0;
	if(ltype == START_BOUNDARY) return 0; 
	if(rtype == END_BOUNDARY) return 0; 

	PSIMI pei = locate_boundary_iterators(*mmap, lpos, rpos);
	SIMI lit = pei.first, rit = pei.second;

	if(lit == mmap->end() || rit == mmap->end()) return 0;

	SIMI it = lit;
	int32_t minc = INT32_MAX;
	int32_t maxc = -1;
	int32_t mins = 0;
	int32_t mint = 0;
	int32_t maxs = 0;
	int32_t maxt = 0;
	while(true)
	{
		if(it->second > maxc)
		{
			maxc = it->second;
			maxs = lower(it->first);
			maxt = upper(it->first);
		}
		if(it->second < minc)
		{
			minc = it->second;
			mins = lower(it->first);
			mint = upper(it->first);
		}

		if(it == rit) break;
		it++;
	}

	//if(rpos - lpos <= 100) return 0;

	printf("large-region %d-%d, len = %d, #intervals = %lu/%lu, minc/maxc = %d/%d, min = %d-%d, max = %d-%d\n", lpos, rpos, rpos - lpos, jmap.size(), boost::distance(jmap.begin(), jmap.end()), minc, maxc, mins, mint, maxs, maxt);
	return 0;
}

int region::smooth_join_interval_map()
{
	vector<PI32> v;
	int32_t p = lpos;
	for(JIMI it = jmap.begin(); it != jmap.end(); it++)
	{
		int32_t p1 = lower(it->first);
		int32_t p2 = upper(it->first);
		assert(p1 >= p);
		assert(p2 > p1);
		if(p1 - p <= cfg.min_subregion_gap) v.push_back(PI32(p, p1));
		p = p2;
	}

	if(p < rpos && rpos - p <= cfg.min_subregion_gap) v.push_back(PI32(p, rpos));

	for(int i = 0; i < v.size(); i++)
	{
		jmap += make_pair(ROI(v[i].first, v[i].second), 1);
	}

	for(JIMI it = jmap.begin(); it != jmap.end(); it++)
	{
		assert(it->second == 1);
	}
	return 0;
}

bool region::empty_subregion(int32_t p1, int32_t p2)
{
	assert(p1 < p2);
	assert(p1 >= lpos && p2 <= rpos);

	//printf(" region = [%d, %d), subregion [%d, %d), length = %d\n", lpos, rpos, p1, p2, p2 - p1);
	if(p2 - p1 < cfg.min_subregion_length) return true;

	PSIMI pei = locate_boundary_iterators(*mmap, p1, p2);
	SIMI it1 = pei.first, it2 = pei.second;
	if(it1 == mmap->end() || it2 == mmap->end()) return true;

	int32_t sum = compute_sum_overlap(*mmap, it1, it2);
	double ratio = sum * 1.0 / (p2 - p1);
	//printf(" region = [%d, %d), subregion [%d, %d), overlap = %.2lf\n", lpos, rpos, p1, p2, ratio);
	//if(ratio < min_subregion_overlap + max_intron_contamination_coverage) return true;
	if(ratio < cfg.min_subregion_overlap) return true;

	return false;
}

int region::build_partial_exons()
{
	pexons.clear();

	//printf("size = %lu, size2 = %lu, [%d, %d), [%d, %d)\n", jmap.size(), distance(jmap.begin(), jmap.end()), lower(jmap.begin()->first), upper(jmap.begin()->first), lpos, rpos);

	assert(lpos < rpos);
	if(jmap.size() == 0 && rpos == lpos + 1 && (ltype == END_BOUNDARY || rtype == START_BOUNDARY))
	{
		partial_exon pe(lpos, rpos, ltype, rtype);
		pe.ave = cfg.min_guaranteed_edge_weight;
		pe.dev = 1.0;
		pexons.push_back(pe);
		return 0;
	}

	if(jmap.size() >= 1 && lower(jmap.begin()->first) == lpos && upper(jmap.begin()->first) == rpos)
	{
		partial_exon pe(lpos, rpos, ltype, rtype);
		evaluate_rectangle(*mmap, pe.lpos, pe.rpos, pe.ave, pe.dev, pe.max);
		pexons.push_back(pe);
		return 0;
	}

	if(ltype == RIGHT_SPLICE && jmap.find(ROI(lpos, lpos + 1)) == jmap.end())
	{
		partial_exon pe(lpos, lpos + 1, ltype, END_BOUNDARY);
		pe.ave = cfg.min_guaranteed_edge_weight;
		pe.dev = 1.0;
		pexons.push_back(pe);
	}

	for(JIMI it = jmap.begin(); it != jmap.end(); it++)
	{
		int32_t p1 = lower(it->first);
		int32_t p2 = upper(it->first);
		assert(p1 < p2);
		
		bool b = empty_subregion(p1, p2);
		if(p1 == lpos && ltype == RIGHT_SPLICE) b = false;
		if(p2 == rpos && rtype == LEFT_SPLICE) b = false;
		if(b == true) continue;

		int lt = (p1 == lpos) ? ltype : START_BOUNDARY;
		int rt = (p2 == rpos) ? rtype : END_BOUNDARY;

		partial_exon pe(p1, p2, lt, rt);
		evaluate_rectangle(*mmap, pe.lpos, pe.rpos, pe.ave, pe.dev, pe.max);
		pexons.push_back(pe);
	}

	if(rtype == LEFT_SPLICE && jmap.find(ROI(rpos - 1, rpos)) == jmap.end())
	{
		partial_exon pe(rpos - 1, rpos, START_BOUNDARY, rtype);
		pe.ave = cfg.min_guaranteed_edge_weight;
		pe.dev = 1.0;
		pexons.push_back(pe);
	}

	return 0;
}

bool region::left_inclusive()
{
	if(pexons.size() == 0) return false;
	if(pexons[0].lpos == lpos) return true;
	else return false;
}

bool region::right_inclusive()
{
	if(pexons.size() == 0) return false;
	if(pexons[pexons.size() - 1].rpos == rpos) return true;
	else return false;
}

int region::print(int index) const
{
	int32_t lc = compute_overlap(*mmap, lpos);
	int32_t rc = compute_overlap(*mmap, rpos - 1);
	printf("region %d: partial-exons = %lu, type = (%d, %d), pos = [%d, %d), boundary coverage = (%d, %d)\n", 
			index, pexons.size(), ltype, rtype, lpos, rpos, lc, rc);

	/*
	for(JIMI it = jmap.begin(); it != jmap.end(); it++)
	{
		printf(" [%d, %d) -> %d\n", lower(it->first), upper(it->first), it->second);
	}
	*/

	return 0;
}

int region::calculate_significance()
{
	if(pexons.size() == 1 && pexons[0].lpos == lpos && pexons[0].rpos == rpos) return 0;

	int read_length = sp.insertsize_median / 2;
	int total_reads = 0;
	int total_length = rpos - lpos;
	for(int i = 0; i < pexons.size(); i++)
	{
		partial_exon &pe = pexons[i];
		int length = pe.rpos - pe.lpos;
		int reads = 1 + pe.ave * length / read_length;
		total_reads += reads;
	}

	vector<double> pvalues(pexons.size(), 1);

	while(true)
	{
		bool flag = false;
		for(int i = 0; i < pexons.size(); i++)
		{
			if(pvalues[i] <= cfg.min_subregion_pvalue) continue;

			partial_exon &pe = pexons[i];
			int length = pe.rpos - pe.lpos;
			int reads = 1 + pe.ave * length / read_length;
			double pr = 1.0 * length / total_length;
			if(pr <= 0) pr = 0;
			if(pr >= 1) pr = 1;

			pvalues[i] = compute_binomial_pvalue(total_reads, pr, reads) * total_length;

			if(cfg.verbose >= 2)
			{
				printf("subregion %d-%d, type = (%d, %d), range = %d-%d, ltype = %d, rtype = %d, total-length = %d, pr = %.4lf, total-reads = %d, reads = %d, pvalue = %.8lf\n",
						pe.lpos, pe.rpos, pe.ltype, pe.rtype, lpos, rpos, ltype, rtype, total_length, pr, total_reads, reads, pvalues[i]);
			}

			if(pvalues[i] > cfg.min_subregion_pvalue) continue;

			flag = true;
			total_length -= length;
			total_reads -= reads;
		}
		if(flag == false) break;
	}

	for(int i = 0; i < pexons.size(); i++)
	{
		partial_exon &pe = pexons[i];
		pe.pvalue = pvalues[i];

		if(pe.ave < cfg.min_subregion_overlap) pe.pvalue = 1;
		if(pe.rpos - pe.lpos < cfg.min_subregion_length) pe.pvalue = 1;
		if(pe.lpos == lpos && ltype == RIGHT_SPLICE) pe.pvalue = 0;
		if(pe.rpos == rpos && rtype == LEFT_SPLICE) pe.pvalue = 0;

		if(cfg.verbose >= 2)
		{
			printf("subregion %d-%d, type = (%d, %d), range = %d-%d, ltype = %d, rtype = %d, pvalue = %.8lf\n", pe.lpos, pe.rpos, pe.ltype, pe.rtype, lpos, rpos, ltype, rtype, pe.pvalue);
		}
	}
	return 0;
}

long double region::calculate_score(int n, int k, int z)
{
	// compute the probability of 
	// x1 + x2 + ... xk = n
	// and there exists one i s.t. xi >= z

	vector<long double> v;
	for(int i = 1; i <= k; i++)
	{
		if(i * z > n) break;
		long double s = log_factorial(n, k, z, i);
		v.push_back(s);
	}

	if(v.size() <= 0) return -9999999;

	long double ans = v[0];
	for(int i = 2; i < v.size(); i += 2) ans = log_add(ans, v[i]);
	for(int i = 1; i < v.size(); i += 2) ans = log_subtract(ans, v[i]);
	printf("ans = %.4Lf, total count = %Lf\n", ans, exp(ans));
	ans -= log_factorial(n, k, 0, k);
	printf("n = %d, k = %d, z = %d, ans = %.3Lf\n", n, k, z, ans);
	return ans;
}

long double region::log_factorial(int n, int k, int z, int i)
{
	// compute (k choose i) (n+k-1-iz choose k-1)
	assert(i >= 1);
	assert(k >= i);
	assert(n >= i * z);
	long double ans = log_select(n - i * z + k - 1, k - 1);
	ans += log_select(k, i);
	printf("n = %d, k = %d, z = %d, i = %d, ans = %.3Lf, count = %.2Lf\n", n, k, z, i, ans, exp(ans));
	return ans;

	/*
	long double ans = log(k);
	for(int j = n - i * z + 1; j <= n - i * z + k - 1; j++)
	{
		ans += log(j);
	}
	for(int j = 1; j <= k - i; j++)
	{
		ans -= log(j);
	}
	for(int j = 1; j <= i; j++)
	{
		ans -= log(j);
	}
	printf("n = %d, k = %d, z = %d, i = %d, ans = %.3Lf, count = %.2Lf\n", n, k, z, i, ans, exp(ans));
	return ans;
	*/
}

long double region::log_add(long double x, long double y)
{
	if(x < y) return log_add(y, x);
	assert(x >= y);
	//printf("compute exp of %.4Lf = %.4Lf\n", y - x, expm1(y - x));
	long f = floor(y);
	double dy = y - f;
	double dx = x - f;
	return x + log1p(1 + expm1(y - x));
	if(x - y >= 1) return x + log1p(1 + expm1(y - x));
	else return x + log1p((1 + expm1(dy)) / (1 + expm1(dx)));
}

long double region::log_subtract(long double x, long double y)
{
	if(x < y)
	{
		printf("x = %.3Lf, y = %.3Lf\n", x, y);
	}
	//printf("compute exp of %.4Lf = %.4Lf\n", y - x, expm1(y - x));

	long f = floor(y);
	double dy = y - f;
	double dx = x - f;
	return x + log1p(-1 - expm1(y - x));
	//if(x - y >= 1) return x + log1p(-1 - expm1(y - x));
	//else return x + log1p(0 - (1 + expm1(dy)) / (1 + expm1(dx)));
}

long double region::log_select(int n, int k)
{
	long double ans = 0;
	for(int j = 1; j <= k; j++)
	{
		ans += log(1.0L * (n - k + j) / j);
	}
	//printf("log-select: %d choose %d log = %.3Lf, count = %.3Lf\n", n, k, ans, exp(ans));
	return ans;
}
