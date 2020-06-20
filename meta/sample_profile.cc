/*
Part of meta-scallop
(c) 2020 by Mingfu Shao, The Pennsylvania State University
See LICENSE for licensing.
*/

#include "sample_profile.h"
#include "htslib/bgzf.h"
#include "constants.h"

mutex sample_profile::bam_lock;
mutex sample_profile::gtf_lock;

sample_profile::sample_profile()
{
	sfn = NULL;
	hdr = NULL;
	bridged_bam = NULL;
	individual_gtf = NULL;
	data_type = DEFAULT;
	insertsize_low = 80;
	insertsize_high = 500;
	insertsize_median = 250;
	library_type = UNSTRANDED;
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
	for(int i = 0; i < hdr->n_targets; i++)
	{
		hts_itr_t *iter = sam_itr_querys(idx, hdr, hdr->target_name[i]);
		iters.push_back(iter);
	}

	hts_idx_destroy(idx);
	close_align_file();
	return 0;
}

int sample_profile::free_index_iterators()
{
	for(int i = 0; i < iters.size(); i++)
	{
		hts_itr_destroy(iters[i]);
	}
	return 0;
}

int sample_profile::print()
{
	printf("file = %s, type = %d\n", align_file.c_str(), data_type);
	return 0;
}
