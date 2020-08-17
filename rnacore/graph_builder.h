/*
Part of aletsch
(c) 2020 by Mingfu Shao, The Pennsylvania State University
See LICENSE for licensing.
*/

#ifndef __GRAPH_BUILDER_H__
#define __GRAPH_BUILDER_H__

#include "junction.h"
#include "region.h"
#include "partial_exon.h"
#include "splice_graph.h"
#include "parameters.h"
#include "bundle_base.h"
#include "sample_profile.h"

using namespace std;

class graph_builder
{
public:
	graph_builder(bundle_base &bd, const parameters &cfg, const sample_profile &sp);

public:
	const parameters &cfg;			// parameters
	const sample_profile &sp;		// sample profile
	bundle_base &bd;				// given bundle_base
	vector<junction> junctions;		// junctions
	vector<region> regions;			// regions
	vector<partial_exon> pexons;	// partial exons
	vector<bool> regional;			// if a partial_exon is regional

public:
	int build(splice_graph &gr);
	int clear();
	int print(int index);

public:
	int analyze_junctions();

private:
	int build_junctions();
	int remove_opposite_junctions();
	int build_regions();
	int build_partial_exons();
	int filter_partial_exons();
	int build_regional();
	int link_partial_exons();
	int build_splice_graph(splice_graph &gr);
};

#endif
