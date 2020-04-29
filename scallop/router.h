/*
Part of Scallop Transcript Assembler
Part of meta-scallop 
(c) 2017 by Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
(c) 2020 by Mingfu Shao, The Pennsylvania State University.
See LICENSE for licensing.
*/

#ifndef __ROUTER_H__
#define __ROUTER_H__

#include <vector>
#include "util.h"
#include "splice_graph.h"
#include "equation.h"
#include "undirected_graph.h"
#include "hyper_set.h"
#include "parameters.h"

typedef pair<int, double> PID;
typedef map<int, double> MID;
typedef pair<PI, double> PPID;
typedef map<PI, double> MPID;

using namespace std;

class router
{
public:
	router(int r, splice_graph &g, MEI &ei, VE &ie, const parameters &cfg);
	router(int r, splice_graph &g, MEI &ei, VE &ie, const MPII &mpi, const parameters &cfg);
	router(int r, splice_graph &g, MEI &ei, VE &ie, const MPII &mpi, const MPII &mpx, const parameters &cfg);
	router& operator=(const router &rt);

public:
	const parameters &cfg;
	int root;					// central vertex
	splice_graph &gr;			// reference splice graph
	MEI &e2i;					// reference map of edge to index
	VE &i2e;					// reference map of index to edge
	vector<PI> routes;			// pairs of connections
	vector<int> counts;			// counts for routes
	vector<PI> routex;			// pairs of connections
	vector<int> countx;			// counts for routes

	MI e2u;						// edge to index
	vector<int> u2e;			// index to edge
	MED u2w;					// weights of edges of ug
	undirected_graph ug;		// bipartite graph

	int type;					// trivial, splitable, single, or multiple 
	int degree;					// level
	double ratio;				// worst ratio
	vector<equation> eqns;		// split results
	MPID pe2w;					// decompose results (for pairs of edges)

public:
	int classify();												// compute status
	int build();												// give solution

	// init
	int build_indices();										// build u2e and e2u
	int build_bipartite_graph();								// build bipartite graph
	vector<double> compute_balanced_weights();					// balanced weights

	// decompose splitable vertex
	int split();												// for splitable vertices

	// decompose unsplitable vertex with greedy algorithm
	int thread();												// for unsplitable vertices
	int thread_isolate1(int k, vector<double> &vw);
	int thread_isolate2(int k, vector<double> &vw);
	bool thread_leaf(vector<double> &vw);
	bool thread_turn(vector<double> &vw);

	PI get_largest_in_route(int x, const vector<PI> &rr, const vector<int> &cc);
	PI get_largest_out_route(int x, const vector<PI> &rr, const vector<int> &cc);

	// print and stats
	int print();
	int stats();
};

bool compare_edge_weight(const PED &x, const PED &y);

#endif
