#ifndef __PREDECOMPOSE_H__
#define __PREDECOMPOSE_H__

#include "splice_graph.h"
#include "hyper_set.h"
#include "equation.h"
#include "router.h"
#include "path.h"

typedef map< edge_descriptor, vector<int> > MEV;
typedef pair< edge_descriptor, vector<int> > PEV;
typedef pair< vector<int>, vector<int> > PVV;
typedef pair<PEE, int> PPEEI;
typedef map<PEE, int> MPEEI;
typedef pair<int, int> PI;
typedef map<int, int> MI;

class predecompose
{
public:
	predecompose(const splice_graph &gr);

public:
	int assemble();
	int write(ostream &os);

public:
	splice_graph gr;					// splice graph
	MEI e2i;							// edge map, from edge to index
	VE i2e;								// edge map, from index to edge
	MEV mev;							// super edges
	vector<int> v2v;					// vertex map
	vector<path> paths;					// predicted paths
	set<int> nonzeroset;				// vertices with degree >= 1

private:
	// init
	int init_vertex_map();
	int init_super_edges();
	int init_inner_weights();
	int init_nonzeroset();

	// decompose
	int resolve_trivial_vertex();
	int decompose_trivial_vertex(int v);
	int decompose_vertex_replace(int v, MPID &pe2w);
	int merge_adjacent_edges(int x, int y);
	int merge_adjacent_edges(int x, int y, double ww);
	int merge_adjacent_equal_edges(int x, int y);
	int split_edge(int exi, double w);
	int remove_edge(int e);
	int collect_paths();
	int collect_path(int e, int s, int t);
	int balance_vertex(int x);
};

#endif
