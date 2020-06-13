/*
Part of meta-scallop
(c) 2020 by Mingfu Shao, The Pennsylvania State University
See LICENSE for licensing.
*/

#ifndef __COMBINED_GROUP_H__
#define __COMBINED_GROUP_H__

#include "parameters.h"
#include "disjoint_set.h"
#include "combined_graph.h"
#include "constants.h"
#include <mutex>

class graph_group
{
public:
	graph_group(string c, char s, const parameters &cfg);

public:
	const parameters &cfg;
	vector<combined_graph> gset;	// given graphs
	vector< vector<int> > gvv;		// merged graphs
	string chrm;
	char strand;

private:
	MISI sindex;				// splice index
	vector<bool> grouped;		// track grouped graphs
	static mutex gmutex;		// global mutex
	double min_similarity;		// minimum similarity for this round
	int min_group_size;			// minimum #graphs to form a group

public:
	int add_graph(const combined_graph &gr);
	int resolve();
	int print();
	int stats(int k);

private:
	int build_splice_index();
	int process_subset1(const set<int> &ss);
	int process_subset2(const set<int> &ss, disjoint_set &ds);
	int build_similarity(const vector<int> &ss, vector<PPID> &vpid, bool local);
	int augment_disjoint_set(const vector<PPID> &vpid, disjoint_set &ds);
	int build_groups(const vector<int> &ss, disjoint_set &ds);
	int build_groups(disjoint_set &ds);
	vector<PPID> filter(const vector<PPID> &vpid);
	vector<PPID> filter(const vector<int> &ss, const vector<PPID> &vpid);
	vector<int> filter(const set<int> &s);
};

#endif
