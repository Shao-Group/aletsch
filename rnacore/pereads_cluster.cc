#include "pereads_cluster.h"
#include "util.h"
#include <cstdio>

pereads_cluster::pereads_cluster()
{
	bounds.assign(5, 0);
	count = 0;
}

int pereads_cluster::clear()
{
	chain1.clear();
	chain2.clear();
	bounds.clear();
	count = 0;
	return 0;
}

int pereads_cluster::print(int index) const
{
	printf("pereads-cluster %d: #reads = %d\n", index, count);
	return 0;
}
