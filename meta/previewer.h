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
	previewer(const parameters &cfg, sample_profile &sp);
	~previewer();

private:
	const parameters &cfg;
	sample_profile &sp;

	int max_preview_reads;
	int max_preview_spliced_reads;
	int min_preview_spliced_reads;
	double preview_infer_ratio;

public:
	int infer_library_type();
	int infer_insertsize();
	int process(bundle &bb, map<int32_t, int> &m);
};

#endif
