#ifndef __ESSENTIAL_H__
#define __ESSENTIAL_H__

#include "splice_graph.h"
#include "hyper_set.h"

using namespace std;

typedef map<vector<int>, int> MVII;

int build_child_splice_graph(splice_graph &root, splice_graph &gr, const set<int> &ss);
int build_child_hyper_set(hyper_set &hyper, hyper_set &hs, const set<int> &ss);
int build_path_coordinates(splice_graph &gr, const vector<int> &v, vector<int32_t> &vv);

#endif
