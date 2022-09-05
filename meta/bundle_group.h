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
#include <mutex>

class bundle_group
{
public:
	bundle_group(string c, char s, const parameters &cfg);
	bundle_group(const bundle_group &p) = default;
	bundle_group(bundle_group &&p) = default;

public:
	const parameters &cfg;
	vector<bundle> gset;				// given graphs
	vector<vector<int32_t>> splices;	// splices for all bundles
	vector<join_interval_map> jmaps;	// join interval maps for all bundles
	vector<vector<int>> gvv;			// merged graphs
	string chrm;
	char strand;

private:
	MISI sindex;				// splice index
	interval_set_map jindex;	// index for jmaps
	vector<bool> grouped;		// track grouped graphs
	static mutex gmutex;		// global mutex
	double min_similarity;		// minimum similarity for this round
	int min_group_size;			// minimum #graphs to form a group

public:
	//int add_graph(const bundle &gr);
	int resolve();
	int print();
	int stats(int k);

private:
	int build_splices();
	int build_join_interval_maps();
	int build_splice_index();
	int build_join_interval_map_index();
	int process_subset1(const set<int> &ss);
	int process_subset2(const set<int> &ss, disjoint_set &ds, int sim);
	int build_splice_similarity(const vector<int> &ss, vector<PPID> &vpid, bool local);
	int build_overlap_similarity(const vector<int> &ss, vector<PPID> &vpid, bool local);
	int augment_disjoint_set(const vector<PPID> &vpid, disjoint_set &ds);
	int build_groups(const vector<int> &ss, disjoint_set &ds);
	int build_groups(disjoint_set &ds);
	int test_overlap_similarity();
	vector<PPID> filter(const vector<PPID> &vpid);
	vector<PPID> filter(const vector<int> &ss, const vector<PPID> &vpid);
	vector<int> filter(const set<int> &s);
};

#endif
