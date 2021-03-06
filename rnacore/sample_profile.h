/*
Part of aletsch
(c) 2020 by Mingfu Shao, The Pennsylvania State University
See LICENSE for licensing.
*/

#ifndef __SAMPLE_PROFILE_H__
#define __SAMPLE_PROFILE_H__

#include <vector>
#include <string>
#include <htslib/sam.h>
#include <mutex>
#include <fstream>

using namespace std;

class sample_profile
{
public:
	sample_profile(int id);

public:
	int sample_id;
	string align_file;
	string index_file;
	samFile *sfn;
	bam_hdr_t *hdr;
	BGZF *bridged_bam;
	ofstream *individual_gtf;
	static mutex bam_lock;
	static mutex gtf_lock;
	int data_type;
	int library_type;
	int bam_with_xs;
	int insertsize_low;
	int insertsize_high;
	int insertsize_median;
	double insertsize_ave;
	double insertsize_std;
	vector<double> insertsize_profile;
	vector<hts_itr_t*> iters;

public:
	int load_profile(const string &dir);
	int save_profile(const string &dir);
	int open_align_file();
	int open_bridged_bam(const string &dir);
	int init_bridged_bam(const string &dir);
	int open_individual_gtf(const string &dir);
	int read_align_headers();
	int read_index_iterators();
	int free_align_headers();
	int free_index_iterators();
	int close_individual_gtf();
	int close_bridged_bam();
	int close_align_file();
	int print();
};

#endif
