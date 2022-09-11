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
#include <thread>
#include <boost/asio/post.hpp>
#include <boost/asio/thread_pool.hpp>
#include <boost/pending/disjoint_sets.hpp>

typedef boost::asio::thread_pool thread_pool;

class assembler
{
public:
	assembler(const parameters &cfg, transcript_set_pool &tspool, mutex &mylock, thread_pool &pool);

public:
	const parameters &cfg;
	thread_pool &pool;
	mutex &mylock;
	transcript_set_pool &tspool;

public:
	int resolve(vector<bundle*> gv, int instance);
	int build_similarity(vector<bundle*> &gv, vector<vector<PID>> &sim);
	int assemble(vector<bundle*> gv, int instance);
	int assemble(bundle &cb, int instance);
	int assemble(splice_graph &gx, phase_set &px, int sid);
	int transform(bundle &cb, splice_graph &gr, bool revising);
	int fix_missing_edges(splice_graph &gr, splice_graph &gx);
	int bridge(vector<bundle*> gv);

	int refine_pairwise(vector<bundle*> gv, vector<vector<PID>> &sim);
	int refine(vector<bundle*> gv);
	int refine(bundle *bd, splice_graph &gr);
	int pairwise_assemble(vector<bundle*> gv, transcript_set &ts, vector<vector<PID>> &sim, int instance);
	int bridge_pairwise(vector<bundle*> gv, vector<vector<PID>> &sim);
};

#endif
