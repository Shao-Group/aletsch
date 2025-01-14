/*
Part of aletsch
(c) 2020 by Mingfu Shao, The Pennsylvania State University
See LICENSE for licensing.
*/

#include "hit.h"
#include "sample_profile.h"
#include "htslib/bgzf.h"
#include "constants.h"
#include "parameters.h"
#include <cassert>
#include <cmath>

mutex sample_profile::bam_lock;
mutex sample_profile::gtf_lock;

sample_profile::sample_profile(int id, int32_t p)
{
	sample_id = id;
	sfn = NULL;
	hdr = NULL;
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

int sample_profile::open_individual_ftr(const string &dir)
{
	char file[10240];
	sprintf(file, "%s/%d.trstFeature.csv", dir.c_str(), sample_id);
	individual_ftr = new ofstream;
	individual_ftr->open(file, std::ofstream::app);
    individual_ftr->setf(ios::fixed, ios::floatfield);
    individual_ftr->precision(2);
	if(individual_ftr->fail()) 
	{
		printf("cannot open individual feature file %s\n", file);
		exit(0);
	}
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

int sample_profile::close_individual_ftr()
{
	individual_ftr->close();
	delete individual_ftr;
	return 0;
}

int sample_profile::close_align_file()
{
    if(hdr != NULL) bam_hdr_destroy(hdr);
    if(sfn != NULL) sam_close(sfn);
	return 0;
}

int sample_profile::set_batch_boundaries(int min_bundle_gap, int max_read_span)
{
	open_align_file();

	start1.resize(hdr->n_targets);
	start2.resize(hdr->n_targets);
	start_off.resize(hdr->n_targets);
	end1.resize(hdr->n_targets);
	end2.resize(hdr->n_targets);

	for(int i = 0; i < hdr->n_targets; i++)
	{
		int32_t len = hdr->target_len[i];
		int n = hdr->target_len[i] / region_partition_length + 1; 
		//printf("hdr size = %d, chrm %d len = %d, n = %d\n", hdr->n_targets, i, len, n);
		start1[i].assign(n, 0);
		start2[i].assign(n, 0);
		start_off[i].assign(n, 0);
		end1[i].assign(n, 0);
		end2[i].assign(n, 0);
	}

	int hid = 0;
	int tid = -1;
	int rid = 0;
	int32_t rpos = 0;
	bam1_t *b1t = bam_init1();
    while(sam_read1(sfn, hdr, b1t) >= 0)
	{
		bam1_core_t &p = b1t->core;

		if((p.flag & 0x4) >= 1) continue;												// read is not mapped
		//if((p.flag & 0x100) >= 1) continue;											// secondary alignment
		//if(p.n_cigar > cfg.max_num_cigar) continue;									// ignore hits with more than max-num-cigar types
		//if(p.qual < cfg.min_mapping_quality) continue;								// ignore hits with small quality
		//if(p.n_cigar < 1) continue;													// should never happen

		hit ht(b1t, hid++);
        if(fabs(ht.pos - ht.rpos) >= max_read_span) continue;
		//ht.set_tags(b1t);
		//ht.set_strand(library_type);

		if(ht.tid != tid)
		{
			if(tid >= 0) end1[tid][rid] = rpos;
			assert(ht.tid < start1.size());
			tid = ht.tid;
			rid = 0;
			off_t offt = bgzf_tell(sfn->fp.bgzf);
			start1[tid][rid] = ht.pos;
			start2[tid][rid] = ht.rpos;
			start_off[tid][rid] = offt;
			rpos = ht.rpos;
		}

		if(ht.pos > rpos + min_bundle_gap)
		{
			if(ht.pos >= region_partition_length * (1 + rid))
			{
				end1[tid][rid] = rpos;
				rid = ht.pos / region_partition_length;
				assert(rid < start1[tid].size());
				off_t offt = bgzf_tell(sfn->fp.bgzf);
				start1[tid][rid] = ht.pos;
				start2[tid][rid] = ht.rpos;
				start_off[tid][rid] = offt;
			}
		}

		if(ht.rpos > rpos) rpos = ht.rpos;
	}

	/*for(int i = 0; i < hdr->n_targets; i++)
	{
		int32_t len = hdr->target_len[i];
		for(int k = 0; k < start1[i].size(); k++)
		{
			printf("boundaries of tid %d, region %d: %d(%d)-%d | %d-%d, len = %d\n", 
					i, k, start1[i][k], start2[i][k], end1[i][k], k * region_partition_length, (k+1)* region_partition_length, len);
		}
	}*/

    bam_destroy1(b1t);
	close_align_file();
	return 0;
}

int sample_profile::read_input_gtf_file()
{
	// TODO: construct input_gtf_trsts and input_gtf_map
	return 0;
}

int sample_profile::read_index_iterators()
{
	open_align_file();
	hts_idx_t *idx = sam_index_load(sfn, index_file.c_str());

	printf("hdr size = %d\n", hdr->n_targets);

	iters.clear();
	iters.resize(hdr->n_targets);
	start1.resize(hdr->n_targets);
	start2.resize(hdr->n_targets);
	end1.resize(hdr->n_targets);
	end2.resize(hdr->n_targets);
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
			printf("build index for target-id %d, %d-%d, query = %s\n", i, s, t, query.c_str());

			//hts_itr_t *iter = sam_itr_querys(idx, hdr, query.c_str());
			hts_itr_t *iter = NULL;
			iters[i].push_back(iter);
			start1[i].push_back(s);
			start2[i].push_back(s);
			end1[i].push_back(s + region_partition_length);
			end2[i].push_back(s + region_partition_length);
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
