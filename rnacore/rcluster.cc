#include "rcluster.h"
#include "util.h"
#include <cstdio>

int rcluster::clear()
{
	vv.clear();
	vl.clear();
	vr.clear();
	return 0;
}

int rcluster::print(int index) const
{
	printf("rcluster %d: #reads = %lu, ", index, vv.size());
	printf("  vv = ( ");
	printv(vv);
	printf(")\n");
	return 0;
}
