/*
Part of aletsch
(c) 2020 by Mingfu Shao, The Pennsylvania State University
See LICENSE for licensing.
*/

#ifndef __INCUBATOR_H__
#define __INCUBATOR_H__

#include "splice_graph.h"
#include "interval_map.h"
#include "bundle.h"
#include "bundle_group.h"
#include "parameters.h"
#include "transcript_set.h"
#include <ctime>
#include <mutex>
#include <thread>
#include <boost/asio/post.hpp>
#include <boost/asio/thread_pool.hpp>
#include <boost/pending/disjoint_sets.hpp>


typedef map< int32_t, set<int> > MISI;
typedef pair< int32_t, set<int> > PISI;
typedef pair<int, int> PI;
typedef boost::asio::thread_pool thread_pool;

class incubator
{
public:
	incubator(vector<parameters> &v);
	~incubator();

public:
	vector<parameters> &params;						// parameters 
	vector<sample_profile> samples;					// samples
	map<string, vector<PI>> sindex;					// sample index
	vector<bundle_group> groups;					// graph groups
	transcript_set_pool tspool;						// a pool for ts
	//vector<transcript_set> tsets;					// transcript sets for instances
	transcript_set tmerge;							// assembled transcripts for all samples
	ofstream meta_gtf;								// meta gtf

public:
	int resolve();

	int generate(string chrm);
	int merge();
	int assemble();
	int rearrange();
	int postprocess();

private:
	int read_bam_list();
	int init_samples();
	int free_samples();
	int build_sample_index();
	int generate(sample_profile &sp, int tid, string chrm, mutex &mylock);
	int assemble(vector<bundle*> gv, int instance, mutex &mylock, thread_pool &pool);
	int postprocess(const transcript_set &ts, ofstream &fout, mutex &mylock);
	int save_transcript_set(const transcript_set &ts, mutex &mylock);
	int write_individual_gtf(int id, const vector<transcript> &vt, const vector<int> &ct, const vector<pair<int, double>> &v);
	int print_groups();
};

#endif
