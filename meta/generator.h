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
	generator(sample_profile &sp, vector<combined_graph> &cbv, const parameters &c);
	~generator();

private:
	const parameters &cfg;
	sample_profile &sp;
	samFile *sfn;
	bam_hdr_t *hdr;
	bam1_t *b1t;

	vector<combined_graph> &vcb;

	int index;
	int qcnt;
	double qlen;

public:
	int resolve();

private:
	int generate(bundle *bb, mutex &mylock, int index);
	int partition(splice_graph &gr, phase_set &hs, vector<pereads_cluster> &ub, vector<splice_graph> &grv, vector<phase_set> &hsv, vector< vector<pereads_cluster> > &ubv);
	bool process_regional_graph(splice_graph &gr, phase_set &ps, vector<pereads_cluster> &vc);
};

#endif
