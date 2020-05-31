/*
Part of meta-scallop
(c) 2020 by Mingfu Shao, The Pennsylvania State University
See LICENSE for licensing.
*/

#ifndef __ASSEMBLER_H__
#define __ASSEMBLER_H__

#include "combined_graph.h"
#include "parameters.h"
#include "sample_profile.h"
#include "transcript_set.h"
#include "splice_graph.h"
#include <mutex>

class assembler
{
public:
	assembler(const parameters &cfg);

public:
	const parameters &cfg;							// parameters 

public:
	int assemble(vector<combined_graph*> gv, int batch, int instance, transcript_set &ts, vector<sample_profile> &samples);
	int resolve_cluster(vector<combined_graph*> gv, combined_graph &cb, vector<sample_profile> &samples);
	int assemble(combined_graph &cb, transcript_set &ts, int mode);
	int assemble(combined_graph &cb, vector<transcript> &vt);
	int assemble(splice_graph &gx, phase_set &px, vector<transcript> &vt);
};

#endif
