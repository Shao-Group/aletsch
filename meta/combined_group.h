#ifndef __COMBINED_GROUP_H__
#define __COMBINED_GROUP_H__

#include "combined_graph.h"
#include <mutex>

typedef map< int32_t, set<int> > MISI;
typedef pair< int32_t, set<int> > PISI;
typedef pair<int, int> PI;
typedef pair<PI, double> PID;
typedef map<int, double> MID;

class combined_group
{
public:
	combined_group(string c, char s);

public:
	vector<combined_graph> gset;		// given graphs
	vector<combined_graph> mset;		// merged graphs
	string chrm;
	char strand;

private:
	MISI mis;
	vector<PID> vpid;

public:
	int add_graph(const combined_graph &gr);
	int resolve(int max_combined, double ratio);
	int write(mutex &mylock, ofstream &fout, int offset);
	int stats();

private:
	int build_splice_map();
	int build_similarity(double ratio);
	int combine_graphs(int max_combined);
};

bool compare_graph_similarity(const PID &x, const PID &y);

#endif
