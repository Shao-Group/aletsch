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
	vector<string> bams;						// bam files
	vector<transcript_set> tss;					// assembled transcripts

	vector<map<string, int>> g2g;				// group map
	vector<combined_group> groups;				// graph groups

public:
	int resolve();
	int read_bam_list();
	int clear();
	int postprocess();

	int generate(int a, int b);
	int merge();
	int assemble();

private:
	int generate(const string &file, vector<combined_group> &gv, mutex &mylock);
	int assemble(vector<combined_graph*> gv, int instance, mutex &mylock);
	int postprocess(const transcript_set &ts, ofstream &fout, mutex &mylock);
	int init_transcript_sets();
	int assemble(vector<combined_graph*> gv, int instance, int subindex, transcript_set &ts);
	int assemble(combined_graph &cb, transcript_set &ts);
	int resolve_cluster(vector<combined_graph*> gv, combined_graph &cb);
	int store_transcripts(const transcript_set &ts, mutex &mylock);

	int print_groups();
};

#endif
