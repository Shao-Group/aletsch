/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#ifndef __PREVIEWER_H__
#define __PREVIEWER_H__

#include "hit.h"
#include "scallop/config.h"
#include "bundle_base.h"
#include "sample_profile.h"

#include <fstream>
#include <string>
#include <map>

using namespace std;

class previewer
{
public:
	previewer(const string &file);
	~previewer();

private:
	string input_file;
	samFile *sfn;
	bam_hdr_t *hdr;
	bam1_t *b1t;

	int max_preview_reads;
	int max_preview_spliced_reads;
	int min_preview_spliced_reads;
	double preview_infer_ratio;

public:
	int open_file();
	int close_file();
	int infer_library_type(config &cfg, sample_profile &sp);
	int infer_insertsize(config &cfg, sample_profile &sp);
	int process(bundle_base &bb, config &cfg, map<int32_t, int> &m);
};

#endif
