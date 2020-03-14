/*
Part of Coral
(c) 2019 by Mingfu Shao and The Pennsylvania State University.
See LICENSE for licensing.
*/

#include "hyper_node.h"
#include <cstdio>

hyper_node::hyper_node(const vector<int> &a, const vector<int> &b, int c)
{
	origin = a;
	todate = b;
	count = c;
}

int hyper_node::print(int index) const
{
	printf("hyper-node %d: count = %d", index, count);
	for(int k = 0; k < origin.size(); k++)
	{
		printf("(%d, %d) ", origin[k], todate[k]);
	}
	return 0;
}
