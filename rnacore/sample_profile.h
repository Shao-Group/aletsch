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
	sample_profile(int id, int32_t p);

public:
	int sample_id;
	string align_file;
	string index_file;
	samFile *sfn;
	bam_hdr_t *hdr;
	ofstream *individual_gtf;
	ofstream *individual_ftr;
	static mutex bam_lock;
	static mutex gtf_lock;
	int data_type;
	int spn;
	int num_xs;
	int library_type;
	int bam_with_xs;
	int insert_total;
	int insertsize_low;
	int insertsize_high;
	int insertsize_median;
	double insertsize_ave;
	double insertsize_std;
	int32_t region_partition_length;
	vector<double> insertsize_profile;
	vector<vector<hts_itr_t*>> iters;
	vector<vector<int32_t>> start1;
	vector<vector<int32_t>> start2;
	vector<vector<int32_t>> end1;
	vector<vector<int32_t>> end2;

public:
	int set_batch_boundaries(int gap);
	int load_profile(const string &dir);
	int save_profile(const string &dir);
	int open_align_file();
	int open_individual_gtf(const string &dir);
	int open_individual_ftr(const string &dir);
	int read_align_headers();
	int read_index_iterators();
	int free_align_headers();
	int free_index_iterators();
	int close_individual_gtf();
	int close_individual_ftr();
	int close_align_file();
	int print();
};

#endif
