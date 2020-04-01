#ifndef __PHASE_SET_H__
#define __PHASE_SET_H__

#include <map>
#include <set>
#include <vector>

#include "util.h"

using namespace std;

typedef pair<vector<int>, int> PVII;
typedef map<vector<int>, int> MVII;

class phase_set
{
public:
	MVII nodes;			// hyper-edges using list-of-nodes

public:
	int clear();
	int add_node_list(const set<int> &s);
	int add_node_list(const set<int> &s, int c);
	int add_node_list(const vector<int> &s, int c);
	int print();
};

#endif
