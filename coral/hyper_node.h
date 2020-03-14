/*
Part of Coral
(c) 2019 by Mingfu Shao and The Pennsylvania State University.
See LICENSE for licensing.
*/

#ifndef __HYPER_NODE_H__
#define __HYPER_NODE_H__

#include <vector>

using namespace std;

class hyper_node
{
public:
	hyper_node(const vector<int> &a, const vector<int> &b, int c);

public:
	vector<int> origin;
	vector<int> todate;
	int count;

public:
	int print(int index) const;
};

#endif
