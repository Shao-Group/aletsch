/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#ifndef __GENERATOR_H__
#define __GENERATOR_H__

#include <fstream>
#include <string>
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
	generator(const string &bamfile, vector<combined_graph> &cbv, parameters &c);
	~generator();

private:
	parameters &cfg;
	sample_profile sp;
	samFile *sfn;
	bam_hdr_t *hdr;
	bam1_t *b1t;
	bundle bb1;		// +
	bundle bb2;		// -
	vector<bundle> pool;

	vector<combined_graph> &vcb;

	int index;
	int qcnt;
	double qlen;

public:
	int resolve();

private:
	int generate(int n);
	int partition(splice_graph &gr, phase_set &hs, const vector<pereads_cluster> &ub, vector<splice_graph> &grv, vector<phase_set> &hsv, vector< vector<pereads_cluster> > &ubv);
};

#endif
