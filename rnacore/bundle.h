/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#ifndef __BUNDLE_H__
#define __BUNDLE_H__

#include <stdint.h>
#include <cstring>
#include <string>
#include <vector>

#include "hit.h"
#include "interval_map.h"

using namespace std;

class bundle
{
public:
	bundle();

public:
	string chrm;					// chromosome name
	char strand;					// strandness
	int32_t tid;					// chromosome ID
	int32_t lpos;					// the leftmost boundary on reference
	int32_t rpos;					// the rightmost boundary on reference
	vector<hit> hits;				// hits
	vector<PI> frgs;				// fragments (pairs of hits)
	vector<bool> brdg;				// if frgs are bridged
	vector<junction> jcns;			// splice junctions
	chain_set hcst;					// chain set for hits 
	chain_set fcst;					// chain set for frgs
	split_interval_map mmap;		// matched interval map
	split_interval_map imap;		// indel interval map

public:
	int add_hit_intervals(const hit &ht, bam1_t *b);
	bool overlap(const hit &ht) const;
	int clear();
	int print(int index);
	int compute_strand(int libtype);
	int check_left_ascending();
	int check_right_ascending();
	int build_fragments();
	int build_junctions();
	int update_bridged_fragments(const vector<int> &frlist, const vector<int32_t> &chain, const vector<int32_t> &whole);

private:
	int add_hit(const hit &ht);
	int add_intervals(bam1_t *b);
};

#endif
