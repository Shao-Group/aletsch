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

using namespace std;

class generator
{
public:
	generator();
	~generator();

private:
	samFile *sfn;
	bam_hdr_t *hdr;
	bam1_t *b1t;
	bundle_base bb1;		// +
	bundle_base bb2;		// -
	vector<bundle_base> pool;

	ofstream grout;
	int index;
	int qcnt;
	double qlen;

public:
	int resolve();

private:
	int process(int n);
	int generate(const splice_graph &gr, const hyper_set &hs);
	int write_graph(splice_graph &gr, hyper_set &hs);
};

#endif
