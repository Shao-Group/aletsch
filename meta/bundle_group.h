/*
Part of aletsch
(c) 2020 by Mingfu Shao, The Pennsylvania State University
See LICENSE for licensing.
*/

#ifndef __BUNDLE_GROUP_H__
#define __BUNDLE_GROUP_H__

#include "parameters.h"
#include "disjoint_set.h"
#include "bundle2.h"
#include "constants.h"
#include <mutex>

class bundle_group
{
public:
	bundle_group(string c, char s, const parameters &cfg);
	bundle_group(const bundle_group &p) = default;
	bundle_group(bundle_group &&p) = default;

public:
	const parameters &cfg;
	vector<bundle2> gset;		// given graphs
	vector<vector<int32_t>> splices;	// splices for all bundles
	vector<vector<int>> gvv;			// merged graphs
	string chrm;
	char strand;

private:
	MISI sindex;				// splice index
	vector<bool> grouped;		// track grouped graphs
	static mutex gmutex;		// global mutex
	double min_similarity;		// minimum similarity for this round
	int min_group_size;			// minimum #graphs to form a group

public:
	//int add_graph(const bundle2 &gr);
	int resolve();
	int print();
	int stats(int k);

private:
	int build_splices();
	int build_splice_index();
	int process_subset1(const set<int> &ss);
	int process_subset2(const set<int> &ss, disjoint_set &ds);
	int build_similarity(const vector<int> &ss, vector<PPID> &vpid, bool local);
	int augment_disjoint_set(const vector<PPID> &vpid, disjoint_set &ds);
	int build_groups(const vector<int> &ss, disjoint_set &ds);
	int build_groups(disjoint_set &ds);
	vector<PPID> filter(const vector<PPID> &vpid);
	vector<PPID> filter(const vector<int> &ss, const vector<PPID> &vpid);
	vector<int> filter(const set<int> &s);
};

#endif
