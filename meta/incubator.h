#ifndef __INCUBATOR_H__
#define __INCUBATOR_H__

#include "splice_graph.h"
#include "interval_map.h"
#include "combined_graph.h"
#include "combined_group.h"
#include "merged_graph.h"
#include "transcript.h"
#include "scallop/config.h"
#include <mutex>

typedef map< int32_t, set<int> > MISI;
typedef pair< int32_t, set<int> > PISI;
typedef pair<int, int> PI;

class incubator
{
public:
	incubator(const config &c);

public:
	config cfg;									// config for scallop
	vector<combined_group> groups;				// graph groups
	vector< map<string, int> > g2g;				// group map
	map< size_t, vector<transcript> > trsts;	// assembled transcripts

public:
	int resolve();

public:
	int load();
	int merge();
	int assemble();
	int postprocess();

	// write and print
	int write(const string &file, bool headers = false);
	int print_groups();
};

int load_multiple(const vector<string> &files, vector<combined_group> &gv, mutex &mylock, vector< map<string, int> > &g2g, const config &cfg);
int load_single(const string &file, vector<combined_graph> &vc);

int assemble_single(combined_graph &cb, int instance, map< size_t, vector<transcript> > &trsts, mutex &mylock, const config &cfg);
int index_transcript(map< size_t, vector<transcript> > &mt, const transcript &t);
bool query_transcript(const map< size_t, vector<transcript> > &mt, const transcript &t);

#endif
