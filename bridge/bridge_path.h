/*
Part of meta-scallop
(c) 2020 by Mingfu Shao, The Pennsylvania State University
See LICENSE for licensing.
*/

#ifndef __bridge_path_H__
#define __bridge_path_H__

#include <vector>
#include <stdint.h>

using namespace std;

class bridge_path
{
public:
	bridge_path();
	~bridge_path();

public:
	int type;
	vector<int> v;
	vector<int32_t> chain;
	vector<int32_t> whole;
	double score;
	vector<int> stack;
	int count;

public:
	int clear();
	int print(int index) const;
	int print_bridge(int index) const;
	vector<int> index(int n) const;
};

bool compare_bridge_path_vertices(const bridge_path &p1, const bridge_path &p2);
bool compare_bridge_path_score(const bridge_path &p1, const bridge_path &p2);
bool compare_bridge_path_stack(const bridge_path &p1, const bridge_path &p2);

#endif
