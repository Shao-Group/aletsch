/*
Part of aletsch
(c) 2020 by Mingfu Shao, The Pennsylvania State University
See LICENSE for licensing.
*/

#ifndef __CHAIN_SET_H__
#define __CHAIN_SET_H__

#include <map>
#include <set>
#include <vector>

#include "util.h"

using namespace std;

class chain_set
{
private:
	map<int, PI> hmap;
	map<int32_t, int> pmap;
	vector<vector<vector<int32_t>>> chains;

public:
	vector<int32_t> get_chain(int h);
	int add(const vector<int32_t> &v, int h);
	int remove(int h);
	int print();
	int clear();
};

#endif
