/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
Part of meta-scallop
(c) 2020 by Mingfu Shao, The Pennsylvania State University
See LICENSE for licensing.
*/

#ifndef __PREVIEWER_H__
#define __PREVIEWER_H__

#include "parameters.h"
#include "hit.h"
#include "bundle.h"
#include "sample_profile.h"

#include <fstream>
#include <string>
#include <map>

using namespace std;

class previewer
{
public:
	previewer(const string &file, const parameters &cfg, sample_profile &sp);
	~previewer();

private:
	string input_file;
	const parameters &cfg;
	sample_profile &sp;
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
	int infer_library_type();
	int infer_insertsize();
	int process(bundle &bb, map<int32_t, int> &m);
};

#endif
