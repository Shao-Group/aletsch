/*
Part of meta-scallop
(c) 2020 by Mingfu Shao, The Pennsylvania State University
See LICENSE for licensing.
*/

#include "sample_profile.h"
#include "htslib/bgzf.h"

sample_profile::sample_profile()
{
	bridged_bam = NULL;
}

int sample_profile::open_bridged_bam(const string &dir)
{
	samFile *sfn = sam_open(file_name.c_str(), "r");
	bam_hdr_t *hdr = sam_hdr_read(sfn);

	char file[10240];
	sprintf(file, "%s/%d.bam", dir.c_str(), sample_id);
	bridged_bam = bgzf_open(file, "w");
	bam_hdr_write(bridged_bam, hdr);

    bam_hdr_destroy(hdr);
    sam_close(sfn);
	return 0;
}

int sample_profile::close_bridged_bam()
{
	if(bridged_bam == NULL) return 0;
	bgzf_close(bridged_bam);
	return 0;
}
