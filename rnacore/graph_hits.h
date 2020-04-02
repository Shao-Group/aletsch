#ifndef __GRAPH_HITS_H__
#define __GRAPH_HITS_H__

#include "splice_graph.h"
#include "pereads_cluster.h"
#include "phase_set.h"
#include "hit.h"

using namespace std;

class graph_hits
{
public:
	graph_hits(splice_graph &gr, vector<hit> &hits);

public:
	splice_graph &gr;				// given splice graph
	vector<hit> &hits;				// given hits

public:
	int group_pereads(vector< vector<PI> > &vc, vector<bool> &paired);
	int build_pereads_clusters(const vector<PI> &fs, vector<pereads_cluster> &vc);
	int build_phase_set_from_unpaired_reads(const vector<bool> &paired, phase_set &ps);
	vector< vector<int> > partition(vector< vector<int32_t> > &fs, int r);
};

bool compare_rank0(const vector<int32_t> &x, const vector<int32_t> &y);
bool compare_rank1(const vector<int32_t> &x, const vector<int32_t> &y);
bool compare_rank2(const vector<int32_t> &x, const vector<int32_t> &y);
bool compare_rank3(const vector<int32_t> &x, const vector<int32_t> &y);

#endif
