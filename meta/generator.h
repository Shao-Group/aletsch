/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#ifndef __GENERATOR_H__
#define __GENERATOR_H__

#include <fstream>
#include <string>
#include "bundle_base.h"
#include "bundle.h"
#include "splice_graph.h"
#include "hyper_set.h"
#include "combined_graph.h"
#include "scallop/config.h"

using namespace std;

class generator
{
public:
	generator(vector<combined_graph> &cbv, const config &c);
	~generator();

private:
	config cfg;
	samFile *sfn;
	bam_hdr_t *hdr;
	bam1_t *b1t;

	ofstream grout;
	vector<combined_graph> &vcb;

	int index;
	int qcnt;
	double qlen;

public:
	int resolve();

private:
	int generate(bundle &bd);
	int write_graph(splice_graph &gr, hyper_set &hs);
};

#endif
