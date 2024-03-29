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

#define INTERVAL_BUF_SIZE 10

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
	vector<int32_t> splices;		// list of splicing positions
	chain_set hcst;					// chain set for hits 
	chain_set fcst;					// chain set for frgs
	split_interval_map mmap;		// matched interval map
	split_interval_map imap;		// indel interval map
	int32_t interval_buf[INTERVAL_BUF_SIZE * 2];
	int32_t interval_cnt[INTERVAL_BUF_SIZE];

public:
	int clear();
	int print(int index);
	int compute_strand(int libtype);
	int check_left_ascending();
	int check_right_ascending();
	int add_hit_intervals(const hit &ht, bam1_t *b);
	int build_fragments();
	int count_unbridged();
	int build_phase_set(phase_set &ps, splice_graph &gr);
	int update_bridges(const vector<int> &frlist, const vector<int32_t> &chain, int strand);
	int filter_multialigned_hits();
	int add_borrowed_path(const vector<int32_t> &p, double w);
	int add_buf_intervals();

private:
	int add_hit(const hit &ht);
	int add_intervals(bam1_t *b);
	int filter_secondary_hits();
	int eliminate_hit(int k);
	int eliminate_bridge(int k);
	bool overlap(const hit &ht) const;
};

#endif
