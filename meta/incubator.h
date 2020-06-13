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
#include "graph_group.h"
#include "parameters.h"
#include "transcript_set.h"
#include <mutex>

typedef map< int32_t, set<int> > MISI;
typedef pair< int32_t, set<int> > PISI;
typedef pair<int, int> PI;

class incubator
{
public:
	incubator(vector<parameters> &v);
	~incubator();

public:
	vector<parameters> &params;						// parameters 
	vector<sample_profile> samples;					// samples
	map<string, vector<PI>> sindex;					// sample index
	vector<graph_group> groups;						// graph groups
	vector<transcript_set> tsets;					// transcript sets for instances
	transcript_set tmerge;							// assembled transcripts for all samples
	ofstream meta_gtf;								// meta gtf

public:
	int resolve();

	int generate(const vector<PI> &v);
	int merge();
	int assemble();
	int rearrange();
	int postprocess();

private:
	int read_bam_list();
	int set_threading();
	int init_samples();
	int free_samples();
	int build_sample_index();
	int generate(sample_profile &sp, int tid, mutex &mylock);
	int assemble(vector<combined_graph*> gv, int instance, mutex &mylock);
	int postprocess(const transcript_set &ts, ofstream &fout, mutex &mylock);
	int move_transcript_set(transcript_set &ts, mutex &mylock);
	int write_individual_gtf(int id, const vector<transcript> &vt, const vector<int> &v);
	int print_groups();
};

#endif
