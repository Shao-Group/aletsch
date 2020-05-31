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
#include <mutex>

class assembler
{
public:
	assembler(vector<sample_profile> &samples, const parameters &cfg);

public:
	vector<sample_profile> &samples;				// samples
	const parameters &cfg;							// parameters 

public:
	int assemble(vector<combined_graph*> gv, int batch, int instance, transcript_set &ts);
	int assemble(combined_graph &cb, transcript_set &ts, int mode);
	int assemble(combined_graph &cb, vector<transcript> &vt);
	int resolve_cluster(vector<combined_graph*> gv, combined_graph &cb);
};

#endif
