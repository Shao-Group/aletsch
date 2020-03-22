#ifndef __COMBINED_GRAPH_H__
#define __COMBINED_GRAPH_H__

#include <vector>
#include <stdint.h>
#include <iostream>
#include <string>

#include "splice_graph.h"
#include "hyper_set.h"
#include "fcluster.h"
#include "interval_map.h"

using namespace std;

typedef pair<int, int> PI;
typedef pair<int32_t, int32_t> PI32;
typedef pair<double, int> DI;
typedef pair<int32_t, DI> PIDI;
typedef pair<PI32, DI> PPDI;
typedef pair< vector<int32_t>, vector<PPDI> > PVDI;
typedef pair< vector<int32_t>, vector<int32_t> > PV32;
typedef map< vector<int32_t>, vector<int32_t> > MV32;

class combined_graph
{
public:
	combined_graph();

public:
	int num_combined;
	string gid;
	string chrm;
	char strand;
	bool parent;

	vector<PPDI> regions;
	vector<PPDI> junctions;
	vector<PIDI> sbounds;
	vector<PIDI> tbounds;
	vector<PV32> phase;
	vector<PV32> reads;
	vector<int32_t> splices;

	vector<combined_graph> children;

	map<int32_t, int> lindex;
	map<int32_t, int> rindex;
	map<int32_t, int32_t> smap;
	map<int32_t, int32_t> tmap;

public:
	// build from gr, hs, and ub
	int build(splice_graph &gr, hyper_set &hs, vector<fcluster> &ub);
	int build_regions(splice_graph &gr);
	int build_start_bounds(splice_graph &gr);
	int build_end_bounds(splice_graph &gr);
	int build_splices_junctions(splice_graph &gr);
	int build_phase(splice_graph &gr, hyper_set &hs);
	int build_reads(splice_graph &gr, vector<fcluster> &ub);

	// combine (only splices)
	int combine(const combined_graph &gt);
	int get_overlapped_splice_positions(const vector<int32_t> &v) const;

	// combine children
	int combine_children();
	int combine_regions(split_interval_map &imap, const combined_graph &gt);
	int combine_junctions(map<PI32, DI> &m, const combined_graph &gt);
	int combine_start_bounds(map<int32_t, DI> &m, const combined_graph &gt);
	int combine_end_bounds(map<int32_t, DI> &m, const combined_graph &gt);
	int combine_phase(MV32 &m, const combined_graph &gt);
	int combine_reads(MV32 &m, const combined_graph &gt);

	// recover splice graph and phasing paths
	int resolve(splice_graph &gr, hyper_set &hs);
	int build_region_index();
	int group_junctions();
	int build_splice_graph(splice_graph &xr);
	int group_start_boundaries(splice_graph &xr);
	int group_end_boundaries(splice_graph &xr);
	int build_phasing_paths(splice_graph &xr, hyper_set &hs);
	PIDI get_leftmost_bound();
	PIDI get_rightmost_bound();

	// get reliable elements
	set<PI32> get_reliable_junctions(int samples, double weight);
	set<int32_t> get_reliable_splices(int samples, double weight);
	set<int32_t> get_reliable_adjacencies(int samples, double weight);
	set<int32_t> get_reliable_start_boundaries(int samples, double weight);
	set<int32_t> get_reliable_end_boundaries(int samples, double weight);

	// mics
	int print(int index);
	int clear();
};

#endif
