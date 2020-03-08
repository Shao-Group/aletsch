/*
Part of Coral
(c) 2019 by Mingfu Shao, The Pennsylvania State University.
See LICENSE for licensing.
*/

#ifndef __HYPER_GRAPH_H__
#define __HYPER_GRAPH_H__

#include <map>
#include <set>
#include <vector>

#include "util.h"
#include "splice_graph.h"

using namespace std;

typedef pair<vector<int>, int> PVII;
typedef map<vector<int>, int> MVII;
typedef map< int, set<int> > MISI;
typedef pair< int, set<int> > PISI;
typedef vector< vector<int> > VVI;
typedef map< pair<int, int>, int> MPII;
typedef pair< pair<int, int>, int> PPII;

class hyper_graph
{
public:
	hyper_graph(const MVII &m);

public:
	vector< vector<int> > nodes;		// list of phasing paths
	vector<int> vcnts;					// counts for each path

private:
	vector< map<int, int> > adjsetx;    // adj map (out)
	vector< map<int, int> > adjsety;    // adj map (in)

public:
	int build_overlap_index();
	PI get_overlap(const vector<int> &vx, const vector<int> &vy);

	int align_paths(const MVII &m);
	int align_path(const vector<int> &x);
	map<int, int> index_path(const vector<int> &x);

	int print_nodes();
	int print_index();

	int keep_maximal_nodes();
	int keep_compatible_nodes(splice_graph &gr);
};

int compare_phasing_paths(const vector<int> &ref, const vector<int> &qry);
bool identical(const vector<int> &x, int x1, int x2, const vector<int> &y, int y1, int y2);

#endif
