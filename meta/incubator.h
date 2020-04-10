#ifndef __INCUBATOR_H__
#define __INCUBATOR_H__

#include "splice_graph.h"
#include "interval_map.h"
#include "combined_graph.h"
#include "combined_group.h"
#include "transcript.h"
#include "parameters.h"
#include <mutex>

typedef map< int32_t, set<int> > MISI;
typedef pair< int32_t, set<int> > PISI;
typedef pair<int, int> PI;

class incubator
{
public:
	incubator(const parameters &c);

public:
	const parameters &cfg;						// parameters for scallop
	vector<combined_group> groups;				// graph groups
	vector< map<string, int> > g2g;				// group map
	map< size_t, vector<transcript> > trsts;	// assembled transcripts

public:
	int resolve();

public:
	int generate();
	int merge();
	int assemble();
	int postprocess();

	int print_groups();
};

int generate_single(const string &file, vector<combined_group> &gv, mutex &mylock, vector< map<string, int> > &g2g, const parameters &cfg);
int assemble_cluster(vector<combined_graph*> gv, int instance, map< size_t, vector<transcript> > &trsts, mutex &mylock, const parameters &cfg);
int assemble_cluster(vector<combined_graph*> gv, int instance, int subindex, vector<transcript> &vt, const parameters &cfg);
int assemble_single(combined_graph &cb, int instance, map< size_t, vector<transcript> > &trsts, mutex &mylock, const parameters &cfg);
int index_transcript(map< size_t, vector<transcript> > &mt, const transcript &t);
bool query_transcript(const map< size_t, vector<transcript> > &mt, const transcript &t);

#endif
