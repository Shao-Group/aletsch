/*
Part of aletsch
(c) 2020 by Mingfu Shao, The Pennsylvania State University
See LICENSE for licensing.
*/

#ifndef __GRAPH_CLUSTER_H__
#define __GRAPH_CLUSTER_H__

#include "splice_graph.h"
#include "pereads_cluster.h"
#include "phase_set.h"
#include "bundle_base.h"

using namespace std;

class graph_cluster
{
public:
	graph_cluster(splice_graph &gr, bundle_base &bd, int max_gap, bool b);

public:
	splice_graph &gr;				// given splice graph
	bundle_base &bd;				// given bundle_base

private:
	vector<vector<int>> groups;
	vector<int32_t> extend;
	int max_partition_gap;
	bool store_hits;

public:
	int build_pereads_clusters(vector<pereads_cluster> &vc);

private:
	int group_pereads();
	int build_pereads_clusters(int g, vector<pereads_cluster> &vc);
	vector<vector<int>> partition(vector<vector<int32_t>> &fs, int r);
};

#endif
