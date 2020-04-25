/*
Part of meta-scallop
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
	graph_cluster(splice_graph &gr, vector<hit> &hits, int max_gap, bool b);

public:
	splice_graph &gr;				// given splice graph
	vector<hit> &hits;				// given hits

private:
	vector<bool> paired;			
	vector<int32_t> extend;
	vector< vector<PI> > groups;
	int max_partition_gap;
	bool store_hits;

public:
	int build_pereads_clusters(vector<pereads_cluster> &vc);
	int build_phase_set_from_unpaired_reads(phase_set &ps);

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
