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

class combined_group
{
public:
	combined_group(string c, char s, const parameters &cfg);

public:
	const parameters &cfg;
	vector<combined_graph> gset;		// given graphs
	vector< vector<int> > gvv;		// merged graphs
	string chrm;
	char strand;

private:
	MISI mis;
	vector<PPID> vpid;

public:
	int add_graph(const combined_graph &gr);
	int resolve();
	int stats();
	int print();

private:
	int build_splice_map();
	int build_similarity();
	int combine_graphs();
};

bool compare_graph_similarity(const PPID &x, const PPID &y);

#endif
