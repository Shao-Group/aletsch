/*
Part of meta-scallop
(c) 2020 by Mingfu Shao, The Pennsylvania State University
See LICENSE for licensing.
*/

#include "sample_profile.h"
#include "htslib/bgzf.h"
#include "constants.h"

mutex sample_profile::bam_lock;

sample_profile::sample_profile()
{
	sfn = NULL;
	hdr = NULL;
	bridged_bam = NULL;
	data_type = DEFAULT;
}

int sample_profile::open_input_file()
{
	sfn = sam_open(file_name.c_str(), "r");
	hdr = sam_hdr_read(sfn);
	return 0;
}

int sample_profile::open_bridged_bam(const string &dir)
{
	char file[10240];
	sprintf(file, "%s/%d.bam", dir.c_str(), sample_id);
	bridged_bam = bgzf_open(file, "w");
	bam_hdr_write(bridged_bam, hdr);
	return 0;
}

int sample_profile::close_bridged_bam()
{
	if(bridged_bam == NULL) return 0;
	int f = bgzf_close(bridged_bam);
	printf("close bridged bam for sample %s with code %d\n", file_name.c_str(), f);
	return 0;
}

int sample_profile::close_input_file()
{
    if(hdr != NULL) bam_hdr_destroy(hdr);
    if(sfn != NULL) sam_close(sfn);
	return 0;
}

int sample_profile::print()
{
	printf("file = %s, type = %d\n", file_name.c_str(), data_type);
	return 0;
}
