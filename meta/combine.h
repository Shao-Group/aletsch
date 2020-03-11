#ifndef __MERGER_H__
#define __MERGER_H__

#include "splice_graph.h"
#include "interval_map.h"
#include "combined_graph.h"
#include "combined_group.h"
#include "merged_graph.h"
#include "transcript.h"
#include <mutex>

typedef map< int32_t, set<int> > MISI;
typedef pair< int32_t, set<int> > PISI;
typedef pair<int, int> PI;

class incubator
{
public:
	incubator(int m, int t);

public:
	vector<combined_group> groups;			// graph groups
	vector< map<string, int> > g2g;

	vector<combined_graph> fixed;			// fixed set of graphs
	int max_combined;						// parameter
	string mdir;							// output dir
	int max_threads;

public:
	// multiple-thread load
	int load(const string &file);
	int merge(double merge_ratio);

	// write and print
	int write(const string &file, bool headers = false);
	int print_groups();
};

int load_multiple(const vector<string> &files, vector<combined_group> &gv, mutex &mylock, vector< map<string, int> > &g2g);
int load_single(const string &file, vector<combined_graph> &vc);

int assemble();
int assemble(merged_graph cm, vector<merged_graph> children, map< size_t, vector<transcript> > &trsts, mutex &mylock);
int index_transcript(map< size_t, vector<transcript> > &mt, const transcript &t);
bool query_transcript(const map< size_t, vector<transcript> > &mt, const transcript &t);

#endif

