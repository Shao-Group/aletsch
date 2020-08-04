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
	int assemble(vector<bundle*> gv, int batch, int instance, transcript_set &ts);
	int resolve_cluster(vector<bundle*> gv, bundle &cb);
	int assemble(bundle &cb, transcript_set &ts, int mode);
	int assemble(bundle &cb, vector<transcript> &vt);
	int assemble(splice_graph &gx, phase_set &px, vector<transcript> &vt, int combined = 1);
};

#endif
