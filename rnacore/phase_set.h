/*
Part of meta-scallop
(c) 2020 by Mingfu Shao, The Pennsylvania State University
See LICENSE for licensing.
*/


#ifndef __PHASE_SET_H__
#define __PHASE_SET_H__

#include <map>
#include <set>
#include <vector>

#include "util.h"

using namespace std;

typedef pair<vector<int32_t>, int> PVII;
typedef map<vector<int32_t>, int> MVII;

class phase_set
{
public:
	phase_set();
	phase_set(const phase_set &pc) = default;
	phase_set(phase_set &&pc) = default;
	phase_set& operator=(phase_set &&pc) = default;
	phase_set& operator=(const phase_set &pc) = default;

public:
	MVII pmap;			// hyper-edges using list-of-coordinates

public:
	int add(const vector<int32_t> &s, int c);
	int combine(const phase_set &ps);
	int project_boundaries(const map<int32_t, int32_t> &smap, const map<int32_t, int32_t> &tmap);
	int project_junctions(const map<PI32, PI32> &jm);
	int clear();
	int print();
};

#endif
