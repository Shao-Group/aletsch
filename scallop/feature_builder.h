/*
Part of aletsch
(c) 2024 by Mingfu Shao, The Pennsylvania State University
See LICENSE for licensing.
*/

#ifndef __FEATURE_BUILDER_H__
#define __FEATURE_BUILDER_H__

#include "path.h"
#include "splice_graph.h"
#include "hyper_set.h"
#include "parameters.h"
#include "essential.h"

using namespace std;

class feature_builder
{
public:
	feature_builder(const parameters &cfg);

public:
	const parameters &cfg;			// parameters

public:
	int build_input_gtf(splice_graph &gr, const vector<transcript> &trsts, const map<int64_t, vector<int>> & tmap, vector<path> &paths);
	int build_transcripts(splice_graph &gr, vector<path> &paths, vector<transcript> &trsts);
    int outputPhasingPath(splice_graph &gr, hyper_set &hs);
	MVII downsamplePhasingPaths(const MVII& paths, size_t targetSize);
    int update_trst_features(splice_graph &gr, transcript &trst, int i, vector<path> &paths);
    int check_junc_relation(const vector<pair<int,int>>& junc1, const vector<pair<int,int>>& junc2);
    int infer_introns(const vector<pair<int, int>>& junc1, const vector<pair<int, int>>& junc2);
    int unique_junc(const vector<path>& paths, int i);
};

#endif
