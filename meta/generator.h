/*
Part of meta-scallop
(c) 2020 by Mingfu Shao, The Pennsylvania State University
See LICENSE for licensing.
*/

#ifndef __GENERATOR_H__
#define __GENERATOR_H__

#include <fstream>
#include <string>
#include <mutex>
#include "bundle.h"
#include "splice_graph.h"
#include "phase_set.h"
#include "pereads_cluster.h"
#include "combined_graph.h"
#include "sample_profile.h"
#include "parameters.h"

using namespace std;

class generator
{
public:
	generator(sample_profile &sp, vector<combined_graph> &cbv, vector<transcript> &trsts, const parameters &c, int target_id);
	~generator();

private:
	const parameters &cfg;
	sample_profile &sp;
	int target_id;

	vector<combined_graph> &vcb;
	vector<transcript> &trsts;

	int index;
	int qcnt;
	double qlen;

public:
	int resolve();

private:
	int generate(bundle &bb, int index);
	int partition(splice_graph &gr, phase_set &hs, vector<pereads_cluster> &ub, vector<splice_graph> &grv, vector<phase_set> &hsv, vector< vector<pereads_cluster> > &ubv);
	bool regional(splice_graph &gr, phase_set &ps, vector<pereads_cluster> &vc);
	bool assemble(splice_graph &gr, phase_set &ps, vector<pereads_cluster> &vc, mutex &tlock);
};

#endif
