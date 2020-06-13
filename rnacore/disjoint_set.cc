#include "disjoint_set.h"

disjoint_set::disjoint_set()
{
	sizes = NULL;
	ranks = NULL;
	parents = NULL;
	ds = NULL;
}

disjoint_set::disjoint_set(int n)
{
	init(n);
}

int disjoint_set::init(int n)
{
	num = n;
	sizes = new int[n];
	ranks = new int[n];
	parents = new int[n];
	ds = new disjoint_set_t(ranks, parents);
	clear();
	return 0;
}

disjoint_set::~disjoint_set()
{
	if(sizes != NULL) delete[] sizes;
	if(ranks != NULL) delete[] ranks;
	if(parents != NULL) delete[] parents;
	if(ds != NULL) delete ds;
}

int disjoint_set::clear()
{
	for(int k = 0; k < num; k++)
	{
		sizes[k] = 1;
		ranks[k] = -1;
		parents[k] = -1;
		ds->make_set(k);
	}
	return 0;
}
