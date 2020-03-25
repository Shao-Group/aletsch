/*
Part of Coral
(c) 2019 by Mingfu Shao, The Pennsylvania State University
See LICENSE for licensing.
*/

#ifndef __BRIDGER_H__
#define __BRIDGER_H__

#include "splice_graph.h"
#include "fragment.h"
#include "fcluster.h"
#include "pier.h"
#include "hyper_set.h"
#include "hit.h"

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
	bridger(splice_graph &gr, vector<hit> &hits);
	bridger(splice_graph &gr, const vector<fcluster> &ub);

public:
	splice_graph &gr;				// given splice graph

	vector<fcluster> fclusters;		// clusters
	vector<pier> piers;				// piers
	vector<bool> bridged;			// whether given hits are bridged
	hyper_set hs;					// constructed hyper-set

	int dp_solution_size;
	int dp_stack_size;
	int32_t length_median;
	int32_t length_low;
	int32_t length_high;

public:
	int resolve();
	int collect_unbridged_fclusters(vector<fcluster> &ub);
	
	int build_fragments(vector<fragment> &fs);
	int build_fclusters(vector<fragment> &fs);
	int build_piers();
	int nominate();
	int vote();
	int build_hyper_set();

	int locate_vertex(int32_t p, int a, int b);
	bool align_hit(const hit &h, vector<int> &vv);
	int dynamic_programming(int k1, int k2, vector< vector<entry> > &table);
	vector<int> update_stack(const vector<int> &v, int s);
	vector< vector<int> > trace_back(int k, const vector< vector<entry> > &table);
	int32_t compute_aligned_length(fragment &fr, const vector<int> &v);
};

#endif
