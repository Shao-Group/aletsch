#ifndef __GRAPH_HITS_H__
#define __GRAPH_HITS_H__

#include "splice_graph.h"
#include "rcluster.h"
#include "hyper_set.h"
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
	int build_paired_reads_clusters(vector<PRC> &prc, vector<bool> &paired);
	int build_hyper_set_from_unpaired_reads(const vector<bool> &paired, hyper_set &hs);
};

#endif
