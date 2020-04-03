/*
Part of Coral
(c) 2019 by Mingfu Shao, The Pennsylvania State University
See LICENSE for licensing.
*/

#ifndef __BRIDGE_SOLVER_H__
#define __BRIDGE_SOLVER_H__

#include "splice_graph.h"
#include "phase_set.h"
#include "pier.h"
#include "pereads_cluster.h"

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

class bridge_solver
{
public:
	bridge_solver(splice_graph &gr, const vector<pereads_cluster> &vc);

public:
	splice_graph &gr;						// given splice graph
	const vector<pereads_cluster> &vc;		// given paired reads clusters

	vector<PI> vpairs;						// vertices for each cluster
	vector<pier> piers;						// piers
	map<PI, int> pindex;					// piers index
	vector<bridge_path> opt;				// optimal bridge path

	int dp_solution_size;
	int dp_stack_size;
	int32_t length_low;
	int32_t length_high;

public:
	int resolve(vector<pereads_cluster> &ub, phase_set &ps);
	int print();

private:
	int build_bridging_vertices();
	int build_piers();
	int build_piers_index();
	int nominate();
	int dynamic_programming(int k1, int k2, vector< vector<entry> > &table);
	vector<int> update_stack(const vector<int> &v, int s);
	vector< vector<int> > trace_back(int k, const vector< vector<entry> > &table);
	int vote();
	int vote(int r, bridge_path &bbp);
	int collect_unbridged_clusters(vector<pereads_cluster> &v);
	int build_phase_set(phase_set &ps);
};

#endif
