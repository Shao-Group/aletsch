/*
Part of Coral package
(c) 2019 by Mingfu Shao, The Pennsylvania State University
See LICENSE for licensing.
*/

#ifndef __CORAL_H__
#define __CORAL_H__

#include <fstream>
#include <htslib/sam.h>
#include <string>
#include "bundle_base.h"
#include "reference.h"

using namespace std;

typedef pair<int, int> PI;

class coral
{
public:
	coral();
	~coral();

private:
	samFile *sfn;
	BGZF *bam_out;
	bam_hdr_t *hdr;
	bam1_t *b1t;
	reference ref;
	bundle_base bb1;			// +
	bundle_base bb2;			// -
	int hid;
	int index;
	vector<bundle_base> pool;	// bundle pool
	map<PI, bam1_t> bmap1;		// bridged reads in + bundle
	map<PI, bam1_t> bmap2;		// bridged reads in - bundle
	map<int, bam1_t> umap1;		// unbridged reads in + bundle
	map<int, bam1_t> umap2;		// unbridged reads in - bundle
	set<string> bset;			// bridged multiple-mapped reads
	vector<bam1_t> uset;		// unbridged reads

public:
	int resolve();
	int process();
	int write_unstrand();
	int write_multiple();
};

#endif
