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
#include "bridge_path.h"
#include "interval_map.h"
#include "constants.h"
#include "parameters.h"

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
	phase_set ps;
	vector<pereads_cluster> vc;

public:
	// copy gid etc from cb
	int copy_meta_information(const combined_graph &cb);
	int set_gid(int instance, int subindex);

	// build from gr, hs, and ub
	int build(splice_graph &gr, const phase_set &ps, const vector<pereads_cluster> &ub);
	int build_regions(splice_graph &gr);
	int build_start_bounds(splice_graph &gr);
	int build_end_bounds(splice_graph &gr);
	int build_splices_junctions(splice_graph &gr);

	// compare combined graphs with splices
	int get_overlapped_splice_positions(const vector<int32_t> &v) const;

	// combine children
	int combine(combined_graph *cb);
	int combine(vector<combined_graph*> &gv);
	int combine_regions(split_interval_double_map &imap) const;
	int combine_junctions(map<PI32, DI> &m) const;
	int combine_start_bounds(map<int32_t, DI> &m) const;
	int combine_end_bounds(map<int32_t, DI> &m) const;

	// append elements to combined graph
	int append(const pereads_cluster &pc, const bridge_path &bbp);
	int append_regions(const pereads_cluster &pc, const bridge_path &bbp);
	int append_junctions(const pereads_cluster &pc, const bridge_path &bbp);

	// recover splice graph and phasing paths
	int build_splice_graph(splice_graph &gr, const parameters &cfg);
	PIDI get_leftmost_bound();
	PIDI get_rightmost_bound();

	// get reliable elements
	set<int32_t> get_reliable_splices(int samples, double weight);

	// mics
	int print(int index);
	int clear();
};

#endif
