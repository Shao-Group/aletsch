/*
Part of aletsch
(c) 2020 by Mingfu Shao, The Pennsylvania State University
See LICENSE for licensing.
*/

#ifndef __GENERATOR_H__
#define __GENERATOR_H__

#include <fstream>
#include <string>
#include <mutex>
#include <boost/asio/post.hpp>
#include <boost/asio/thread_pool.hpp>
#include <boost/pending/disjoint_sets.hpp>
#include "bundle.h"
#include "splice_graph.h"
#include "phase_set.h"
#include "pereads_cluster.h"
#include "bundle.h"
#include "sample_profile.h"
#include "transcript_set.h"
#include "parameters.h"

using namespace std;

typedef boost::asio::thread_pool thread_pool;

class generator
{
public:
	generator(sample_profile &sp, vector<bundle> &cbv, const parameters &c, thread_pool &pool, int target_id, int region_id);
	~generator();

private:
	const parameters &cfg;
	thread_pool &pool;
	sample_profile &sp;
	int target_id;
	int region_id;
	mutex vcb_mutex;
	vector<bundle> &vcb;
	int index;

public:
	int resolve();

private:
	int process();
	int generate(bundle_base &bb, int index);
	int partition(splice_graph &gr, phase_set &hs, vector<pereads_cluster> &ub, vector<splice_graph> &grv, vector<phase_set> &hsv, vector< vector<pereads_cluster> > &ubv);
	bool regional(splice_graph &gr, phase_set &ps, vector<pereads_cluster> &vc);
	//bool assemble_single(splice_graph &gr, phase_set &ps, vector<pereads_cluster> &vc);
};

#endif
