#ifndef __ESSENTIAL_H__
#define __ESSENTIAL_H__

#include "hit.h"
#include "splice_graph.h"
#include "hyper_set.h"

using namespace std;

typedef map<vector<int>, int> MVII;

vector<int> project_vector(const vector<int> &v, const map<int, int> &m);
int transform_vertex_set_map(const set<int> &s, map<int, int> &m);
int build_child_splice_graph(splice_graph &root, splice_graph &gr, map<int, int> &ss);
int build_child_hyper_set(hyper_set &hyper, hyper_set &hs, map<int, int> &ss);
int build_path_coordinates(splice_graph &gr, const vector<int> &v, vector<int32_t> &vv);
bool check_continue_vertices(splice_graph &gr, int x, int y);
bool check_valid_path(splice_graph &gr, const vector<int> &vv);
bool align_hit_to_splice_graph(const hit &h, splice_graph &gr, vector<int> &vv);
int build_paired_fragments(const vector<hit> &hits, vector<PI> &fs);

#endif
