/*
Part of aletsch
(c) 2020 by Mingfu Shao, The Pennsylvania State University
See LICENSE for licensing.
*/

#ifndef __INCUBATOR_H__
#define __INCUBATOR_H__

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
	map<pair<string, char>, transcript_set> tts;	// transcripts for each chrm
	ofstream meta_gtf;								// meta gtf
	ofstream meta_ftr;								// meta feature
	vector<bundle_group> grps;						// bundle groups
	thread_pool tpool;
	int group_size;									// number of regions in a group
	vector<mutex> gmutex;							// mutex for writing to gset in each bundle_group
	vector<mutex> tmutex;							// mutex for transcripts in each bundle_group
	//transcript_set_pool tspool;					// a pool for ts
	//transcript_set tmerge;						// assembled transcripts for all samples
	//mutex tlock;									// global lock for transcripts

public:
	int resolve();

	int merge();
	int assemble();
	//int rearrange();
	int postprocess();

private:
	int read_bam_list();
	int init_samples();
	int free_samples();
	int build_sample_index();
	int init_bundle_groups();
	int init_transcript_sets();
	int get_max_region(string chrm);
	int get_chrm_index(string chrm, int sid);
	set<int> get_target_list(int sid);
	int get_bundle_group(string chrm, int gid);
	int generate_merge_assemble(string chrm, int gid);
	int generate(int sid, int tid, int rid, string chrm, mutex &curlock);
	int assemble(bundle_group &g, int gid, int gi);
	int write_individual_gtf(int id, const vector<transcript> &t);
	int write_individual_gtf(int sid);
	int write_combined_gtf();
	int print_groups(const vector<bundle_group> &grps);
	//int postprocess(const transcript_set &ts, ofstream &fout, mutex &mylock);
	//int save_transcript_set(const transcript_set &ts, mutex &mylock);
};

#endif
