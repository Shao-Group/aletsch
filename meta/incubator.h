#ifndef __INCUBATOR_H__
#define __INCUBATOR_H__

#include "splice_graph.h"
#include "interval_map.h"
#include "combined_graph.h"
#include "combined_group.h"
#include "parameters.h"
#include "transcript_set.h"
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
	transcript_set tset;						// assembled transcripts

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

int assemble_cluster(vector<combined_graph*> gv, int instance, transcript_set &ts, mutex &mylock, const parameters &cfg);
int assemble_single(combined_graph &cb, int instance, transcript_set &ts, mutex &mylock, const parameters &cfg);
int assemble_cluster(vector<combined_graph*> gv, int instance, int subindex, transcript_set &ts, const parameters &cfg);
int assemble_single(combined_graph &cb, transcript_set &ts, const parameters &cfg, bool group_boundary);
int resolve_cluster(vector<combined_graph*> gv, combined_graph &cb, const parameters &cfg);
int resolve_cluster(vector<combined_graph*> gv, const parameters &cfg);

#endif
