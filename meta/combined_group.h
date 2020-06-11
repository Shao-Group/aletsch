/*
Part of meta-scallop
(c) 2020 by Mingfu Shao, The Pennsylvania State University
See LICENSE for licensing.
*/

#ifndef __COMBINED_GROUP_H__
#define __COMBINED_GROUP_H__

#include "parameters.h"
#include "combined_graph.h"
#include "constants.h"
#include <mutex>
#include "boost/pending/disjoint_sets.hpp"

typedef disjoint_sets<int*, int*> dss;

class combined_group
{
public:
	combined_group(string c, char s, const parameters &cfg);

public:
	const parameters &cfg;
	vector<combined_graph> gset;	// given graphs
	vector< vector<int> > gvv;		// merged graphs
	string chrm;
	char strand;

private:
	int divisors;
	MISI sindex;				// splicindex
	vector<int*> sizes;			// size for dss
	vector<int*> ranks;			// rank for dss
	vector<int*> parents;		// parent for dss
	vector<dss*> dsss;			// list of dss
	vector<bool> bmap;			// if a graph is clustered
	static mutex gmutex;		// global mutex

public:
	int add_graph(const combined_graph &gr);
	int resolve();
	int stats();
	int print();

private:
	int init_disjoint_sets();
	int build_splice_index();
	int build_disjoint_sets();
	int process_subset(const vector<int> &ss);
	int build_similarity(const vector<int> &ss, vector<PPID> &vpid);
	int build_clusters(const vector<int> &ss, const vector<PPID> &vpid);

};

#endif
