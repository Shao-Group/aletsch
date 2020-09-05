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
#include "chain_set.h"
#include "phase_set.h"
#include "splice_graph.h"

using namespace std;

class bundle_base
{
public:
	bundle_base();

public:
	string chrm;					// chromosome name
	char strand;					// strandness
	int32_t tid;					// chromosome ID
	int32_t lpos;					// the leftmost boundary on reference
	int32_t rpos;					// the rightmost boundary on reference
	vector<hit> hits;				// hits
	vector<AI3> frgs;				// fragments <hit1, hit2, type>, type: -1: cannot be bridged; 0: to-be-bridged; 1: bridge with empty; 2: bridge with extra splices
	vector<AI6> grps;				// <bounds, count, type> of grouped paired reads
	chain_set hcst;					// chain set for hits 
	chain_set fcst;					// chain set for frgs
	chain_set pcst;					// chain set for first end of grps
	chain_set qcst;					// chain set for second end of grps
	chain_set gcst;					// chain set for bridges of grps
	split_interval_map mmap;		// matched interval map
	split_interval_map imap;		// indel interval map

public:
	int clear();
	int print(int index);
	int compute_strand(int libtype);
	int check_left_ascending();
	int check_right_ascending();
	int add_hit_intervals(const hit &ht, bam1_t *b);
	int build_fragments();
	int build_phase_set(phase_set &ps, splice_graph &gr);
	int update_bridges(const vector<int> &frlist, const vector<int32_t> &chain);
	int filter_multialigned_hits();

private:
	int add_hit(const hit &ht);
	int add_intervals(bam1_t *b);
	int filter_secondary_hits();
	int eliminate_hit(int k);
	int eliminate_bridge(int k);
	bool overlap(const hit &ht) const;
};

#endif
