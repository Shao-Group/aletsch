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
	incubator(const vector<parameters> &v);
	~incubator();

public:
	const vector<parameters> &params;				// parameters 
	vector<sample_profile> samples;					// samples
	map<string, vector<PI>> sindex;					// sample index
	vector<transcript_set> tsets;					// transcript sets for instances
	vector<transcript_set> tsave;					// assembled transcripts for all samples
	vector<vector<transcript>> strsts;				// predicted transcripts for each sample
	vector<combined_group> groups;					// graph groups

public:
	int resolve();

	int generate(const vector<PI> &v);
	int merge();
	int assemble();
	int rearrange();
	int postprocess();
	int write();

	int read_bam_list();
	int init_samples();
	int build_sample_index();
	int close_samples();
	int write(int id);

	int print_groups();

private:
	int init_sample(sample_profile &sp);
	int generate(sample_profile &sp, int tid, mutex &mylock);
	int assemble(vector<combined_graph*> gv, int instance, mutex &mylock);
	int rearrange(transcript_set &root, const vector<int> &v);
	int postprocess(const transcript_set &ts, ofstream &fout, mutex &mylock);
	int save_transcripts(const vector<transcript> &v, int sid, mutex &mylock);
	int save_transcript_set(const transcript_set &ts, mutex &mylock);
	int move_transcript_set(transcript_set &ts, mutex &mylock);
};

#endif
