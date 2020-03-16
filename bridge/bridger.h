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

public:
	splice_graph &gr;				// given splice graph
	vector<hit> &hits;				// given hits

public:
	map<int32_t, int> lindex;		// index for left-vertex of gr
	map<int32_t, int> rindex;		// index for right-vertex of gr

	vector<fragment> fragments;		// fragments needs to be bridged
	vector<fcluster> fclusters;		// clusters
	vector<pier> piers;				// piers

	int32_t length_median;
	int32_t length_low;
	int32_t length_high;

public:
	int resolve();
	
private:
	int build_vertex_index();
	int build_fragments();
	int build_fclusters();
	int build_piers();
	int bridge();
	int vote();

	int locate_vertex(int32_t p, int a, int b);
	bool align_hit(const hit &h, vector<int> &vv);
	int dynamic_programming(int k1, int k2, vector< vector<entry> > &table);
	vector<int> update_stack(const vector<int> &v, int s);
	vector< vector<int> > trace_back(int k, const vector< vector<entry> > &table);
	int32_t compute_aligned_length(fragment &fr, const vector<int> &v);
};

#endif
