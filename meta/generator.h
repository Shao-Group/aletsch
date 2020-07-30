/*
Part of aletsch
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
#include "transcript_set.h"
#include "parameters.h"

using namespace std;

class generator
{
public:
	generator(sample_profile &sp, vector<combined_graph> &cbv, transcript_set &ts, const parameters &c, int target_id);
	~generator();

private:
	const parameters &cfg;
	sample_profile &sp;
	int target_id;

	vector<combined_graph> &vcb;
	transcript_set &ts;
	int index;

public:
	int resolve();

private:
	int generate(bundle &bb, int index);
	int bridge(bundle &bb);
	int partition(splice_graph &gr, phase_set &hs, vector<pereads_cluster> &ub, vector<splice_graph> &grv, vector<phase_set> &hsv, vector< vector<pereads_cluster> > &ubv);
	int process_large(vector<pereads_cluster> &vc);
	bool regional(splice_graph &gr, phase_set &ps, vector<pereads_cluster> &vc);
	bool assemble_single(splice_graph &gr, phase_set &ps, vector<pereads_cluster> &vc);
	bool assemble_large(splice_graph &gr, phase_set &ps, vector<pereads_cluster> &vc);
};

#endif
