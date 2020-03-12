#ifndef __GRAPH_REVISE_H__
#define __GRAPH_REVISE_H__

#include "config.h"
#include "splice_graph.h"

using namespace std;

typedef pair<int32_t, int32_t> PI32;

int revise_splice_graph(splice_graph &gr, config *cfg);
int revise_splice_graph_full(splice_graph &gr, config *cfg);

VE compute_maximal_edges(splice_graph &gr);
bool extend_boundaries(splice_graph &gr);
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

#endif
