/*
Part of meta-scallop
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
	sample_profile();

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
	int library_type;
	int data_type;
	double insertsize_ave;
	double insertsize_std;
	int insertsize_low;
	int insertsize_high;
	int insertsize_median;
	vector<double> insertsize_profile;
	vector<hts_itr_t*> iters;

public:
	int open_align_file();
	int open_bridged_bam(const string &dir);
	int open_individual_gtf(const string &dir);
	int close_individual_gtf();
	int close_bridged_bam();
	int close_align_file();
	int build_index_iterators();
	int destroy_index_iterators();
	int print();
};

#endif
