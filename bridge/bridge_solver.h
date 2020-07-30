/*
Part of aletsch
(c) 2020 by Mingfu Shao, The Pennsylvania State University
See LICENSE for licensing.
*/

#ifndef __BRIDGE_SOLVER_H__
#define __BRIDGE_SOLVER_H__

#include "splice_graph.h"
#include "phase_set.h"
#include "pier.h"
#include "pereads_cluster.h"
#include "parameters.h"

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
	bridge_solver(splice_graph &gr, vector<pereads_cluster> &vc, const parameters &cfg, int low, int high);

public:
	splice_graph &gr;						// given splice graph
	vector<pereads_cluster> &vc;			// given paired reads clusters
	vector<bridge_path> opt;				// optimal bridge path
	const parameters &cfg;

private:
	vector<edge_descriptor> adjedges;		// added edges
	vector<PI> vpairs;						// vertices for each cluster
	vector<pier> piers;						// piers
	vector<int> bounds;						// groups of piers
	map<PI, int> pindex;					// piers index

	int32_t length_low;
	int32_t length_high;

public:
	int print();
	int build_phase_set(phase_set &ps);
	int collect_unbridged_clusters(vector<pereads_cluster> &v);
	int write_bridged_reads(BGZF *fout);

private:
	int add_adjacent_edges();
	int remove_adjacent_edges();
	int build_bridging_vertices();
	bool check_left_relaxing(const pereads_cluster &pc, int v);
	bool check_right_relaxing(const pereads_cluster &pc, int v);
	int build_piers();
	int build_piers_index();
	int build_bounds();
	int nominate();
	int nominate(int strand);
	int refine_pier(pier &p);
	int dynamic_programming(int k1, int k2, vector< vector<entry> > &table, int strand);
	vector<int> update_stack(const vector<int> &v, int s);
	vector< vector<int> > trace_back(int k, const vector< vector<entry> > &table);
	int vote();
	int vote(int r, bridge_path &bbp);
};

int add_phases_from_bridged_pereads_cluster(const pereads_cluster &pc, const bridge_path &bbp, phase_set &ps);
int add_phases_from_unbridged_pereads_cluster(const pereads_cluster &pc, phase_set &ps);

#endif
