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
	map<int, AI3> hmap;
	map<int32_t, int> pmap;
	vector<vector<PVI3>> chains;

public:
	int add(const vector<int32_t> &v, int h, int xs);	// if h < 0, don't store the handle
	int remove(int h);									// remove handle and decrease count
	int clear();										// clear everything
	int print();										// print
	vector<int32_t> get_chain(int h) const;				// get chain
	vector<int32_t> get_splices() const;				// get the set of all splices
	PVI3 get(int h) const;								// get chain and return count
};

#endif
