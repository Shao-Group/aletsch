/*
Part of meta-scallop
(c) 2020 by Mingfu Shao, The Pennsylvania State University
See LICENSE for licensing.
*/

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
	~incubator();

public:
	const parameters &cfg;							// parameters for scallop
	vector<sample_profile> samples;					// samples
	vector<vector<transcript>> strsts;				// predicted transcripts for each sample
	vector<transcript_set> tsets;					// assembled transcripts for all samples
	vector<map<string, int>> g2g;					// group map
	vector<combined_group> groups;					// graph groups

public:
	int resolve();
	int read_bam_list();
	int clear();
	int postprocess();
	int write();

	int generate(int a, int b);
	int merge();
	int assemble();
	int write(int id);

private:
	int generate(sample_profile &sp, vector<combined_group> &gv, mutex &mylock);
	int assemble(vector<combined_graph*> gv, int instance, mutex &mylock);
	int postprocess(const transcript_set &ts, ofstream &fout, mutex &mylock);
	int init_transcript_sets();
	int assemble(combined_graph &cb, transcript_set &ts, int mode);
	int resolve_cluster(vector<combined_graph*> gv, combined_graph &cb, mutex &mylock);
	int store_transcripts(const transcript_set &ts, mutex &mylock);

	int print_groups();
};

#endif
