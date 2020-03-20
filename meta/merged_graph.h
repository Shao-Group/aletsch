#ifndef __COMBINED_GRAPH_H__
#define __COMBINED_GRAPH_H__

#include <vector>
#include <stdint.h>
#include <iostream>
#include <string>

#include "combined_graph.h"
#include "splice_graph.h"
#include "hyper_set.h"
#include "config.h"

using namespace std;

typedef pair<int, int> PI;
typedef pair<int32_t, int32_t> PI32;
typedef pair<double, int> DI;
typedef pair<int32_t, DI> PIDI;
typedef pair<PI32, DI> PPDI;
typedef pair< vector<int32_t>, vector<PPDI> > PVDI;

class merged_graph
{
public:
	merged_graph();

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
	vector<PVDI> phase;

	map<int32_t, int> lindex;
	map<int32_t, int> rindex;

	map<int32_t, int32_t> smap;
	map<int32_t, int32_t> tmap;

	splice_graph gr;
	hyper_set hs;

public:
	int solve();
	int print(int index);
	int build(combined_graph &gr, string s);
	set<PI32> get_reliable_junctions(int samples, double weight);
	set<int32_t> get_reliable_splices(int samples, double weight);
	set<int32_t> get_reliable_adjacencies(int samples, double weight);
	set<int32_t> get_reliable_start_boundaries(int samples, double weight);
	set<int32_t> get_reliable_end_boundaries(int samples, double weight);
	int clear();
	
private:
	int build_region_index();
	int group_junctions();
	int build_splice_graph(splice_graph &xr);
	int group_start_boundaries(splice_graph &xr);
	int group_end_boundaries(splice_graph &xr);
	int group_phasing_paths();
	int build_phasing_paths();
	bool continue_vertices(int x, int y, splice_graph &xr);
	PIDI get_leftmost_bound();
	PIDI get_rightmost_bound();
	int locate_left_region(int32_t p, int kl, int kr);
	int locate_right_region(int32_t p, int kl, int kr);
};

#endif
