/*
Part of aletsch
(c) 2020 by Mingfu Shao, The Pennsylvania State University
See LICENSE for licensing.
*/

#ifndef __ASSEMBLER_H__
#define __ASSEMBLER_H__

#include "bundle.h"
#include "parameters.h"
#include "transcript_set.h"
#include "splice_graph.h"
#include <mutex>

class assembler
{
public:
	assembler(const parameters &cfg);

public:
	const parameters &cfg;

public:
	int resolve(vector<bundle*> gv, transcript_set &ts, int instance);
	int assemble(bundle &cb, transcript_set &ts, int instance);
	int assemble(vector<bundle*> gv, transcript_set &ts, int instance);
	int assemble(splice_graph &gx, phase_set &px, transcript_set &ts, int sid);
	int refine_pairwise(bundle &cx, bundle &cy);
	int refine_pairwise(splice_graph &gx, splice_graph &gr);
	int transform(bundle &cb, splice_graph &gr, bool revising);
	int fix_missing_edges(splice_graph &gr, splice_graph &gx);
	int bridge(vector<bundle*> gv);
};

#endif
