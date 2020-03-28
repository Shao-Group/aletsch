#ifndef __ESSENTIAL_H__
#define __ESSENTIAL_H__

#include "hit.h"
#include "splice_graph.h"
#include "hyper_set.h"
#include "rcluster.h"

using namespace std;

typedef map<vector<int>, int> MVII;

// map integers
vector<int> project_vector(const vector<int> &v, const map<int, int> &m);

// transform a set to a map
int transform_vertex_set_map(const set<int> &s, map<int, int> &m);

// determine whether a->b is continuous in graph
bool check_continuous_vertices(splice_graph &gr, int x, int y);

// determine whether vv is a valid path in graph
bool check_valid_path(splice_graph &gr, const vector<int> &vv);

// construct sub-splice-graph
int build_child_splice_graph(splice_graph &root, splice_graph &gr, map<int, int> &ss);

// transform between paths and coordinates
int build_exon_coordinates_from_path(splice_graph &gr, const vector<int> &v, vector<int32_t> &vv);
bool build_path_from_exon_coordinates(splice_graph &gr, const vector<int32_t> &v, vector<int> &vv);
bool build_path_from_intron_coordinates(splice_graph &gr, const vector<int32_t> &v, vector<int> &vv);
bool build_path_from_mixed_coordinates(splice_graph &gr, const vector<int32_t> &v, vector<int> &vv);

// align hit to splice graph
bool align_hit_to_splice_graph(const hit &h, splice_graph &gr, vector<int> &vv);

// transform to paths from coordinates
bool transform_to_paths(splice_graph &gr, PRC &prc);

// match paired-end reads
int build_paired_reads(const vector<hit> &hits, vector<PI> &fs);

#endif
