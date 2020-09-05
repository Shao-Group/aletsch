/*
Part of aletsch
(c) 2020 by Mingfu Shao, The Pennsylvania State University
See LICENSE for licensing.
*/

#ifndef __GRAPH_CLUSTER_H__
#define __GRAPH_CLUSTER_H__

#include "splice_graph.h"
#include "phase_set.h"
#include "bundle_base.h"

using namespace std;

class graph_cluster
{
public:
	graph_cluster(splice_graph &gr, bundle_base &bd, int max_gap);

public:
	splice_graph &gr;				// given splice graph
	bundle_base &bd;				// given bundle_base

private:
	vector<vector<int>> groups;
	int max_partition_gap;

private:
	int group_pereads();
	int build_clusters(int g);
	vector<vector<int>> partition(vector<vector<int32_t>> &fs, int r);
};

bool compare_rank0(const vector<int32_t> &x, const vector<int32_t> &y);
bool compare_rank1(const vector<int32_t> &x, const vector<int32_t> &y);
bool compare_rank2(const vector<int32_t> &x, const vector<int32_t> &y);
bool compare_rank3(const vector<int32_t> &x, const vector<int32_t> &y);

#endif
