/*
Part of aletsch
(c) 2020 by Mingfu Shao, The Pennsylvania State University
See LICENSE for licensing.
*/

#ifndef __BUNDLE_GROUP_H__
#define __BUNDLE_GROUP_H__

#include "parameters.h"
#include "disjoint_set.h"
#include "bundle.h"
#include "constants.h"
#include "interval_map.h"
#include "transcript_set.h"
#include <mutex>
#include <unordered_map>
#include <unordered_set>
#include <boost/asio/post.hpp>
#include <boost/asio/thread_pool.hpp>
#include <boost/pending/disjoint_sets.hpp>

typedef boost::asio::thread_pool thread_pool;

class bundle_group
{
public:
	bundle_group(string c, char s, int r, const parameters &cfg, thread_pool &pool);
	//bundle_group(const bundle_group &p) = default;
	//bundle_group(bundle_group &&p) = default;

public:
	const parameters &cfg;
	vector<bundle> gset;				// given graphs
	vector<vector<int32_t>> splices;	// splices for all bundles
	vector<join_interval_map> jmaps;	// join interval maps for all bundles
	vector<vector<int>> gvv;			// merged graphs
	string chrm;						// chrm name
	char strand;						// strandness
	int rid;							// group id
	thread_pool &tpool;					// thread pool
	transcript_set_pool tspool;			// transcripts for this bundle_group
	transcript_set tmerge;				// merged transcripts

private:
	MISI sindex;				// splice index
	interval_set_map jindex;	// index for jmaps
	vector<bool> grouped;		// track grouped graphs

public:
	int resolve();
	int print();
	int clear();
	int stats(int k);

private:
	int build_splices();
	int build_join_interval_maps();
	int build_splice_index();
	int build_join_interval_map_index();
	int build_splice_similarity(const vector<int> &ss, vector<PPID> &vpid, bool local, double d);
	int build_overlap_similarity(const vector<int> &ss, vector<PPID> &vpid, bool local);
	int test_overlap_similarity();
	int process_subset(const set<int> &ss, disjoint_set &ds, double d);
	int augment_disjoint_set(const vector<PPID> &vpid, disjoint_set &ds);
	int build_groups(disjoint_set &ds);
	int filter(const set<int> &s, disjoint_set &ds, vector<int> &v);
};

#endif
