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
	vector<vector<PVI>> chains;

public:
	int add(const vector<int32_t> &v, int h);	// if h < 0, don't store the handle
	int remove(int h);							// remove handle and decrease count
	PVI get(int h) const;						// get chain and return count
	int clear();								// clear everything
	int print();								// print
};

#endif
