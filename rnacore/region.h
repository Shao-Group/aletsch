/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#ifndef __REGION_H__
#define __REGION_H__

#include <stdint.h>
#include <vector>
#include "interval_map.h"
#include "partial_exon.h"
#include "parameters.h"
#include "sample_profile.h"

using namespace std;

class region
{
public:
	region(int32_t _lpos, int32_t _rpos, int _ltype, int _rtype, const split_interval_map *_mmap, const split_interval_map *_imap, const parameters &cfg, const sample_profile &sp);
	~region();

public:
	const parameters &cfg;			// config
	const sample_profile &sp;		// sample profile
	int32_t lpos;					// the leftmost boundary on reference
	int32_t rpos;					// the rightmost boundary on reference
	int ltype;						// type of the left boundary
	int rtype;						// type of the right boundary
	const split_interval_map *mmap;	// pointer to match interval map
	const split_interval_map *imap;	// pointer to indel interval map
	join_interval_map jmap;			// subregion intervals
	vector<partial_exon> pexons;	// generated partial exons

public:
	int print(int index) const;
	bool left_inclusive();
	bool right_inclusive();
	static long double calculate_score(int n, int k, int z);

private:
	bool empty_subregion(int32_t p1, int32_t p2);
	int build_join_interval_map();
	int smooth_join_interval_map();
	int build_partial_exons();
    int build_indel_coverage();
	int calculate_significance();
	static long double log_factorial(int n, int k, int z, int i);
	static long double log_add(long double x, long double y);
	static long double log_subtract(long double x, long double y);
	static long double log_select(int n, int k);
};

#endif
