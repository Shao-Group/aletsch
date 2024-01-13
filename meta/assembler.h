/*
Part of aletsch
(c) 2020 by Mingfu Shao, The Pennsylvania State University
See LICENSE for licensing.
*/

#ifndef __ASSEMBLER_H__
#define __ASSEMBLER_H__

#include "bundle.h"
#include "parameters.h"
#include "transcript_set.h"
#include "splice_graph.h"
#include <mutex>
#include <unordered_map>
#include <boost/asio/post.hpp>
#include <boost/asio/thread_pool.hpp>
#include <boost/pending/disjoint_sets.hpp>

typedef boost::asio::thread_pool thread_pool;

class assembler
{
public:
	assembler(const parameters &cfg, transcript_set &tmerge, mutex &mylock, thread_pool &p, int rid, int gid, int instance);

public:
	const parameters &cfg;
	thread_pool &pool;
	transcript_set &tmerge;
	mutex &mylock;
	int rid;
	int gid;
	int instance;

public:
	int resolve(vector<bundle*> gv);
	int build_similarity(vector<bundle*> &gv, vector<vector<PID>> &sim);
	int assemble(vector<bundle*> gv);
	int assemble(bundle &cb);
	int assemble(splice_graph &gx, phase_set &px, int sid);
	int transform(bundle &cb, splice_graph &gr, bool revising);
	int fix_missing_edges(splice_graph &gr, splice_graph &gx);
	int bridge(vector<bundle*> gv);
	int combine_bundles(bundle &bd, vector<bundle*> gv);

    //sample support
    //int junction_support(int sample_id, splice_graph &gr, splice_graph &gx);
    int junction_support(splice_graph &gr, unordered_map<int64_t, set<int> > &junc2sup, unordered_map<int64_t, unordered_map<int, double>> &sup2abd);
    int start_end_support(int sample_id, splice_graph &gr, splice_graph &gx);
	int start_end_support(vector<splice_graph*> &grv, const vector<int> &idv);
    int non_splicing_support(int sample_id, splice_graph &gr, splice_graph &gx);
    int boundary_extend(int sample_id, splice_graph &gr, splice_graph &gx, int pos_type);

};

#endif
