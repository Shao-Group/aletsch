/*
Part of Coral
(c) 2019 by Mingfu Shao, The Pennsylvania State University
See LICENSE for licensing.
*/

#ifndef __BRIDGER_H__
#define __BRIDGER_H__

#include "splice_graph.h"
#include "rcluster.h"
#include "pier.h"
#include "hyper_set.h"

using namespace std;

class entry
{
public:
	vector<int> stack;
	int32_t length;
	int trace1;
	int trace2;

public:
	int print();
};

bool entry_compare(const entry &x, const entry &y);

class bridger
{
public:
	bridger(splice_graph &gr, const vector<PRC> &vpr);

public:
	splice_graph &gr;				// given splice graph
	const vector<PRC> &vpr;			// given paired reads clusters

	vector<pier> piers;				// piers
	map<PI, int> pindex;			// piers index
	vector<phase> opt;				// optimal bridging paths

	int dp_solution_size;
	int dp_stack_size;
	int32_t length_median;
	int32_t length_low;
	int32_t length_high;

public:
	int resolve();
	int build_piers();
	int build_piers_index();
	int nominate();
	int dynamic_programming(int k1, int k2, vector< vector<entry> > &table);
	vector<int> update_stack(const vector<int> &v, int s);
	vector< vector<int> > trace_back(int k, const vector< vector<entry> > &table);
	int vote();
	int vote(const PRC &prc, phase &bbp);
	int collect_unbridged_clusters(vector<PRC> &v);
	int build_hyper_set(hyper_set &hs);
};

#endif
