/*
Part of Coral
(c) 2019 by Mingfu Shao, The Pennsylvania State University.
See LICENSE for licensing.
*/

#ifndef __BUNDLE2_H__
#define __BUNDLE2_H__

#include "bundle_base.h"
#include "junction.h"
#include "region.h"
#include "fragment.h"
#include "transcript.h"
#include "hyper_node.h"

using namespace std;

class bundle2 : public bundle_base
{
public:
	bundle2(const bundle_base &bb);
	virtual ~bundle2();

public:
	vector<junction> junctions;			// splice junctions
	vector<region> regions;				// regions
	vector<fragment> fragments;			// fragments (identical ones are merged)
	vector<transcript> ref_trsts;		// overlaped genes in reference
	vector< vector<int> > ref_phase;	// phasing paths for ref transcripts
	vector< vector<PI> > ref_index;		// the set of trsts that contain each region
	vector<hyper_node> hpnodes;			// hyper nodes

public:
	virtual int build();
	int print(int index);

public:
	int build_junctions();
	int extend_junctions();
	int build_regions();
	int link_regions();

	int align_hits_transcripts();
	int align_hit(const map<int32_t, int> &m, const hit &h, vector<int> &v);
	int align_transcript(const map<int32_t, int> &m, const transcript &t, vector<int> &v);
	int index_references();
	int locate_region(int32_t x);

	int build_fragments();
	int group_fragments();
	int bridge_fragments();
	int32_t compute_aligned_length(int32_t k1l, int32_t k2r, const vector<int>& v);

	int build_hyper_nodes();

	// write to bam
	int build_bam1_t(bam1_t &b1t, const vector<int> &v, const hit &h, int32_t p1, int32_t p2, char c);
	int add_cigar_match(bam1_t &b1t, int32_t p1, int32_t p2);
	int add_cigar_skip(bam1_t &b1t, int32_t p1, int32_t p2);
	int write_bam(BGZF *fout, map<PI, bam1_t> &bmap, map<int, bam1_t> &umap);
};

#endif
