#include "pereads_cluster.h"
#include "util.h"
#include <cstdio>

pereads_cluster::pereads_cluster()
{
	bounds.assign(4, 0);
	extend.assign(4, 0);
	count = 0;
}

int pereads_cluster::clear()
{
	chain1.clear();
	chain2.clear();
	bounds.clear();
	extend.clear();
	count = 0;
	return 0;
}

int pereads_cluster::print(int index) const
{
	printf("pereads-cluster %d: #reads = %d, bounds = ( ", index, count);
	printv(bounds);
	printf("), extend = ( ");
	printv(extend);
	printf("), chain1 = ( ");
	printv(chain1);
	printf("), chain2 = ( ");
	printv(chain2);
	printf(")\n");
	return 0;
}
