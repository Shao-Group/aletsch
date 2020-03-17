/*
Part of Coral, an efficient tool to bridge mate pairs
(c) 2018 by Mingfu Shao and The Pennsylvania State University.
See LICENSE for licensing.
*/

#include "fcluster.h"
#include "util.h"
#include <cstdio>

int fcluster::clear()
{
	frset.clear();
	v1.clear();
	v2.clear();
	return 0;
}

int fcluster::print(int index) const
{
	printf("fcluster %d: #fragments = %lu, ", index, frset.size());
	printf("  v1 = ( ");
	printv(v1);
	printf("), v2 = ( ");
	printv(v2);
	printf(")\n");
	return 0;
}

bool compare_fcluster_v1_v2(const fcluster &fx, const fcluster &fy)
{
	for(int k = 0; k < fx.v1.size() && k < fy.v1.size(); k++)
	{
		if(fx.v1[k] < fy.v1[k]) return true;
		if(fx.v1[k] > fy.v1[k]) return false;
	}

	if(fx.v1.size() < fy.v1.size()) return true;
	if(fx.v1.size() > fy.v1.size()) return false;

	// NOTE: reverse sorting for v2 (not for v1)
	for(int k = 0; k < fx.v2.size() && k < fy.v2.size(); k++)
	{
		if(fx.v2[k] > fy.v2[k]) return true;
		if(fx.v2[k] < fy.v2[k]) return false;
	}

	if(fx.v2.size() > fy.v2.size()) return true;
	if(fx.v2.size() < fy.v2.size()) return false;

	return false;
}
