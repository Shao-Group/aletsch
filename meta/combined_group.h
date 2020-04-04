#ifndef __COMBINED_GROUP_H__
#define __COMBINED_GROUP_H__

#include "combined_graph.h"
#include <mutex>

class combined_group
{
public:
	combined_group(string c, char s);

public:
	vector<combined_graph> gset;		// given graphs
	vector< set<int> > clusters;		// merged graphs
	string chrm;
	char strand;

private:
	MISI mis;
	vector<PPID> vpid;

public:
	int add_graph(const combined_graph &gr);
	int resolve();
	int stats();

private:
	int build_splice_map();
	int build_similarity();
	int combine_graphs();
};

bool compare_graph_similarity(const PPID &x, const PPID &y);

#endif
