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
#include "hit.h"

using namespace std;

class graph_cluster
{
public:
	graph_cluster(const splice_graph &gr, const vector<hit> &hits, const vector<PI> &frags, int max_gap, bool b);

public:
	const splice_graph &gr;						// given splice graph
	const vector<hit> &hits;					// given hits
	const vector<PI> &frags;					// given paired hits

private:
	vector<bool> paired;			
	vector<int32_t> extend;
	vector< vector<PI> > groups;
	int max_partition_gap;
	bool store_hits;

public:
	int build_pereads_clusters(vector<pereads_cluster> &vc);
	vector<bool> get_paired();

private:
	int group_pereads();
	int build_pereads_clusters(int g, vector<pereads_cluster> &vc);
	vector< vector<int> > partition(vector< vector<int32_t> > &fs, int r);
};

bool compare_rank0(const vector<int32_t> &x, const vector<int32_t> &y);
bool compare_rank1(const vector<int32_t> &x, const vector<int32_t> &y);
bool compare_rank2(const vector<int32_t> &x, const vector<int32_t> &y);
bool compare_rank3(const vector<int32_t> &x, const vector<int32_t> &y);

#endif
