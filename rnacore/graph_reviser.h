/*
Part of aletsch
(c) 2020 by Mingfu Shao, The Pennsylvania State University
See LICENSE for licensing.
*/


#ifndef __GRAPH_REVISE_H__
#define __GRAPH_REVISE_H__

#include "parameters.h"
#include "splice_graph.h"
#include "bundle_base.h"

using namespace std;

typedef pair<int32_t, int32_t> PI32;

int revise_splice_graph(splice_graph &gr, const parameters &cfg);
int revise_splice_graph_full(splice_graph &gr, const parameters &cfg);

VE compute_maximal_edges(splice_graph &gr);
bool extend_boundaries(splice_graph &gr);
bool remove_trivial_vertices(splice_graph &gr);
bool remove_small_junctions(splice_graph &gr);
bool remove_small_exons(splice_graph &gr, int min_exon);
bool remove_inner_boundaries(splice_graph &gr);
bool remove_intron_contamination(splice_graph &gr, double ratio);
bool keep_surviving_edges(splice_graph &gr, double surviving);
int keep_surviving_edges(splice_graph &gr, const set<PI32> &js, double surviving);
int keep_surviving_edges(splice_graph &gr, const set<int32_t> &ps, double surviving);
int keep_surviving_edges(splice_graph &gr, const set<int32_t> &ps, const set<int32_t> &aj, double surviving);
int keep_surviving_edges(splice_graph &gr, const set<int32_t> &ps, const set<int32_t> &aj, const set<int32_t> &sb, const set<int32_t> &tb, double surviving);

int filter_start_boundaries(splice_graph &gr, const set<int32_t> &ps, double surviving);
int filter_end_boundaries(splice_graph &gr, const set<int32_t> &ps, double surviving);
int filter_junctions(splice_graph &gr, const set<int32_t> &ps, double surviving);
int filter_graph(splice_graph &gr, const set<int32_t> &ps, const set<int32_t> &aj, const set<int32_t> &sb, const set<int32_t> &tb, double surviving);

int refine_splice_graph(splice_graph &gr);

int group_start_boundaries(splice_graph &gr, map<int32_t, int32_t> &smap, int32_t max_group_boundary_distance);
int group_end_boundaries(splice_graph &gr, map<int32_t, int32_t> &tmap, int32_t max_group_boundary_distance);

int identify_boundaries(splice_graph &gr, const parameters &cfg);
bool identify_start_boundary(splice_graph &gr, double min_ratio);
bool identify_end_boundary(splice_graph &gr, double min_ratio);
int left_continuous_extend(splice_graph &gr, int x);
int right_continuous_extend(splice_graph &gr, int x);
int add_distant_in_vertices(splice_graph &gr, int x, set<int> &s);
int add_distant_out_vertices(splice_graph &gr, int x, set<int> &s);
int determine_start_boundary(splice_graph &gr, int a, int b, double &max, double &sum);
int determine_end_boundary(splice_graph &gr, int a, int b, double &max, double &sum);

int remove_false_boundaries(splice_graph &gr, bundle_base &bb, const parameters &cfg);
int catch_false_boundaries(splice_graph &gr, bundle_base &bb, const parameters &cfg);

#endif
