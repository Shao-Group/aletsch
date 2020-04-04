#ifndef __COMBINED_GRAPH_H__
#define __COMBINED_GRAPH_H__

#include <vector>
#include <stdint.h>
#include <iostream>
#include <string>

#include "sample_profile.h"
#include "splice_graph.h"
#include "phase_set.h"
#include "pereads_cluster.h"
#include "interval_map.h"
#include "constants.h"

using namespace std;

class combined_graph
{
public:
	combined_graph();

public:
	int num_combined;
	string gid;
	string chrm;
	char strand;
	sample_profile sp;

	vector<PPDI> regions;
	vector<PPDI> junctions;
	vector<PIDI> sbounds;
	vector<PIDI> tbounds;
	vector<int32_t> splices;
	phase_set phases;
	vector<pereads_cluster> preads;

public:
	// build from gr, hs, and ub
	int build(splice_graph &gr, const phase_set &ps, const vector<pereads_cluster> &ub);
	int build_regions(splice_graph &gr);
	int build_start_bounds(splice_graph &gr);
	int build_end_bounds(splice_graph &gr);
	int build_splices_junctions(splice_graph &gr);

	//int combine_extra_bridged_reads(const vector< vector<int32_t> > &exon_chains, const vector<int> &weights);
	//int combine(const combined_graph &gt);

	int get_overlapped_splice_positions(const vector<int32_t> &v) const;

	// combine children
	int combine();
	int combine_regions(split_interval_double_map &imap) const;
	int combine_junctions(map<PI32, DI> &m) const;
	int combine_start_bounds(map<int32_t, DI> &m) const;
	int combine_end_bounds(map<int32_t, DI> &m) const;

	// recover splice graph and phasing paths
	int build_splice_graph(splice_graph &gr);
	PIDI get_leftmost_bound();
	PIDI get_rightmost_bound();

	// get reliable elements
	set<int32_t> get_reliable_splices(int samples, double weight);

	// mics
	int print(int index);
	int clear();
};

#endif
