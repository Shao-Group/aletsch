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
	MVII pmap;			// hyper-edges using list-of-nodes

public:
	int add(const vector<int32_t> &s, int c);
	int clear();
};

#endif
