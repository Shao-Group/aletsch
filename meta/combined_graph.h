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

using namespace std;

typedef pair<int, int> PI;
typedef pair<int32_t, int32_t> PI32;
typedef pair<double, int> DI;
typedef pair<int32_t, DI> PIDI;
typedef pair<PI32, DI> PPDI;

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

	vector<combined_graph> children;

public:
	// build from gr, hs, and ub
	int build(splice_graph &gr, const phase_set &ps, const vector<pereads_cluster> &ub);
	int build_regions(splice_graph &gr);
	int build_start_bounds(splice_graph &gr);
	int build_end_bounds(splice_graph &gr);
	int build_splices_junctions(splice_graph &gr);

	// add extra reads
	//int combine_extra_bridged_reads(const vector< vector<int32_t> > &exon_chains, const vector<int> &weights);

	// combine (only splices)
	int combine(const combined_graph &gt);
	int get_overlapped_splice_positions(const vector<int32_t> &v) const;

	// combine children
	int combine_children();
	int combine_regions(split_interval_map &imap, const combined_graph &gt);
	int combine_junctions(map<PI32, DI> &m, const combined_graph &gt);
	int combine_start_bounds(map<int32_t, DI> &m, const combined_graph &gt);
	int combine_end_bounds(map<int32_t, DI> &m, const combined_graph &gt);

	// recover splice graph and phasing paths
	int build_splice_graph(splice_graph &gr);
	int build_phase_set(phase_set &ps);
	PIDI get_leftmost_bound();
	PIDI get_rightmost_bound();

	// get reliable elements
	set<int32_t> get_reliable_splices(int samples, double weight);

	// mics
	int print(int index);
	int clear();
};

#endif
