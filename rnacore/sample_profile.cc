/*
Part of aletsch
(c) 2020 by Mingfu Shao, The Pennsylvania State University
See LICENSE for licensing.
*/

#include "sample_profile.h"
#include "htslib/bgzf.h"
#include "constants.h"
#include <cassert>

mutex sample_profile::bam_lock;
mutex sample_profile::gtf_lock;

sample_profile::sample_profile(int id, int32_t p)
{
	sample_id = id;
	sfn = NULL;
	hdr = NULL;
	bridged_bam = NULL;
	individual_gtf = NULL;
	data_type = DEFAULT;
	insertsize_low = 80;
	insertsize_high = 500;
	insertsize_median = 250;
	region_partition_length = p;
	library_type = UNSTRANDED;
	bam_with_xs = 0;
	num_xs = 0;
	spn = 0;
	insert_total = 0;
}

int sample_profile::load_profile(const string &dir)
{
	char file[10240];
	sprintf(file, "%s/%d.profile", dir.c_str(), sample_id);
	ifstream fin(file);
	if(fin.fail())
	{
		printf("cannot open profile to read: %s\n", file);
		return 0;
	}

	char line[10240];
	while(fin.getline(line, 10240, '\n'))
	{
		if(string(line) == "") continue;
		stringstream sstr(line);
		char key[10240];
		sstr >> key;
		
		if(string(key) == "library_type") sstr >> library_type;
		if(string(key) == "bam_with_xs") sstr >> bam_with_xs;
		if(string(key) == "insertsize_low") sstr >> insertsize_low;
		if(string(key) == "insertsize_high") sstr >> insertsize_high;
		if(string(key) == "insertsize_median") sstr >> insertsize_median;
		if(string(key) == "insertsize_ave") sstr >> insertsize_ave;
		if(string(key) == "insertsize_std") sstr >> insertsize_std;
	}

	fin.close();
	return 0;
}

int sample_profile::save_profile(const string &dir)
{
	char file[10240];
	sprintf(file, "%s/%d.profile", dir.c_str(), sample_id);
	ofstream fout(file);
	if(fout.fail())
	{
		printf("cannot open profile to write: %s\n", file);
		return 0;
	}

	fout << "library_type" << " " << library_type << endl;
	fout << "bam_with_xs" << " " << bam_with_xs << endl;

	if(data_type == PAIRED_END)
	{
		fout << "insertsize_low" << " " << insertsize_low << endl;
		fout << "insertsize_high" << " " << insertsize_high << endl;
		fout << "insertsize_median" << " " << insertsize_median << endl;
		fout << "insertsize_ave" << " " << insertsize_ave << endl;
		fout << "insertsize_std" << " " << insertsize_std << endl;
	}

	fout.close();
	return 0;
}

int sample_profile::read_align_headers()
{
	open_align_file();
	hdr = sam_hdr_read(sfn);
	close_align_file();
	return 0;
}

int sample_profile::free_align_headers()
{
    if(hdr != NULL) bam_hdr_destroy(hdr);
	return 0;
}

int sample_profile::open_align_file()
{
	sfn = sam_open(align_file.c_str(), "r");
	hdr = sam_hdr_read(sfn);
	return 0;
}

int sample_profile::init_bridged_bam(const string &dir)
{
	char file[10240];
	sprintf(file, "%s/%d.bam", dir.c_str(), sample_id);
	open_align_file();
	bridged_bam = bgzf_open(file, "w");
	bam_hdr_write(bridged_bam, hdr);
	close_bridged_bam();
	close_align_file();
	return 0;
}

int sample_profile::open_bridged_bam(const string &dir)
{
	char file[10240];
	sprintf(file, "%s/%d.bam", dir.c_str(), sample_id);
	bridged_bam = bgzf_open(file, "a");
	return 0;
}

int sample_profile::open_individual_gtf(const string &dir)
{
	char file[10240];
	sprintf(file, "%s/%d.gtf", dir.c_str(), sample_id);
	individual_gtf = new ofstream;
	individual_gtf->open(file, std::ofstream::app);
	if(individual_gtf->fail()) 
	{
		printf("cannot open individual gtf %s\n", file);
		exit(0);
	}
	return 0;
}

int sample_profile::close_individual_gtf()
{
	individual_gtf->close();
	delete individual_gtf;
	return 0;
}

int sample_profile::close_bridged_bam()
{
	if(bridged_bam != NULL) bgzf_close(bridged_bam);
	return 0;
}

int sample_profile::close_align_file()
{
    if(hdr != NULL) bam_hdr_destroy(hdr);
    if(sfn != NULL) sam_close(sfn);
	return 0;
}

int sample_profile::read_index_iterators()
{
	open_align_file();
	hts_idx_t *idx = sam_index_load(sfn, index_file.c_str());

	iters.clear();
	iters.resize(hdr->n_targets);
	start1.resize(hdr->n_targets);
	start2.resize(hdr->n_targets);
	for(int i = 0; i < hdr->n_targets; i++)
	{
		int32_t len = hdr->target_len[i];
		int n = hdr->target_len[i] / region_partition_length + 1; 
		for(int k = 0; k < n; k++)
		{
			int32_t s = (k + 0) * region_partition_length;
			int32_t t = (k + 2) * region_partition_length;
			//string query = string(hdr->target_name[i]) + ":" + to_string(s) + "-" + to_string(t);
			string query = string(hdr->target_name[i]) + ":" + to_string(s);

			//printf("build index for target-id %d, %d-%d, query = %s\n", i, s, t, query.c_str());

			hts_itr_t *iter = sam_itr_querys(idx, hdr, query.c_str());

			iters[i].push_back(iter);
			start1[i].push_back(s);
			start2[i].push_back(s);
		}
		assert(iters[i].size() == start1[i].size());
		assert(iters[i].size() == start2[i].size());
	}

	hts_idx_destroy(idx);
	close_align_file();
	return 0;
}

int sample_profile::free_index_iterators()
{
	for(int i = 0; i < iters.size(); i++)
	{
		for(int k = 0; k < iters[i].size(); k++)
		{
			hts_itr_destroy(iters[i][k]);
		}
	}
	return 0;
}

int sample_profile::print()
{
	printf("file = %s, type = %d\n", align_file.c_str(), data_type);
	return 0;
}
