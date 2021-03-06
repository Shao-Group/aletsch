/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#ifndef __HYPER_SET_H__
#define __HYPER_SET_H__

#include <map>
#include <set>
#include <vector>

#include "util.h"
#include "splice_graph.h"
#include "phase_set.h"

using namespace std;

typedef pair<vector<int>, int> PVII;
typedef map<vector<int>, int> MVII;
typedef map< int, set<int> > MISI;
typedef pair< int, set<int> > PISI;
typedef vector< vector<int> > VVI;
typedef map< pair<int, int>, int> MPII;
typedef pair< pair<int, int>, int> PPII;

class hyper_set
{
public:
	hyper_set();
	hyper_set(splice_graph &gr, const phase_set &ps);

public:
	MVII nodes;			// hyper-edges using list-of-nodes
	VVI edges;			// hyper-edges using list-of-edges
	vector<int> ecnts;	// counts for edges
	MISI e2s;			// index: from edge to hyper-edges

public:
	int clear();
	int add_node_list(const vector<int> &s, int c, int o = 1);
	int build(directed_graph &gr, MEI &e2i);
	int build_edges(directed_graph &gr, MEI &e2i);
	int build_index();
	int update_index();
	set<int> get_intersection(const vector<int> &v);
	MI get_successors(int e);
	MI get_predecessors(int e);
	MPII get_routes(int x, directed_graph &gr, MEI &e2i);
	int print_nodes();
	int print_edges();
	int write(ostream &os) const;

	int filter();
	int filter_nodes(splice_graph &gr);
	int compare(const hyper_set &hx);
	int merge(const hyper_set &hx);
	int merge_node_list(const vector<int> &s, int c);
	int extend(const hyper_set &hx);

	int build_confident_nodes(int min_count);

public:
	int replace(int x, int e);
	int replace(int x, int y, int e);
	int replace(int x, int y, int x2, int y2);
	int replace(const vector<int> &x, int e);
	int replace_strange(const vector<int> &x, int e);
	int remove(int e);
	int remove(const vector<int> &x);
	int remove(const set<int> &x);
	int remove_pair(int x, int y);
	int insert_between(int x, int y, int e);
	int left_break(int e);
	int right_break(int e);
	bool extend(int e);
	bool left_extend(int e);
	bool left_extend(const vector<int> &s);
	bool right_extend(int e);
	bool right_extend(const vector<int> &s);
	bool left_dominate(int e);
	bool right_dominate(int e);
	bool useful(const vector<int> &v, int k1, int k2);
};

#endif
