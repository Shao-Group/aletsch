/*
Part of meta-scallop
(c) 2020 by Mingfu Shao, The Pennsylvania State University
See LICENSE for licensing.
*/


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

int pereads_cluster::project_junctions(const map<PI32, PI32> &jm)
{
	project_chain(jm, chain1);
	project_chain(jm, chain2);

	if(chain1.size() >= 2 && bounds[0] >= chain1.front()) bounds[0] = chain1.front() - 1;
	if(chain1.size() >= 2 && bounds[1] <= chain1.back()) bounds[1] = chain1.back() + 1;
	if(chain2.size() >= 2 && bounds[2] >= chain2.front()) bounds[2] = chain2.front() - 1;
	if(chain2.size() >= 2 && bounds[3] <= chain2.back()) bounds[3] = chain2.back() + 1;

	if(chain1.size() >= 2 && extend[0] >= chain1.front()) extend[0] = chain1.front() - 1;
	if(chain1.size() >= 2 && extend[1] <= chain1.back()) extend[1] = chain1.back() + 1;
	if(chain2.size() >= 2 && extend[2] >= chain2.front()) extend[2] = chain2.front() - 1;
	if(chain2.size() >= 2 && extend[3] <= chain2.back()) extend[3] = chain2.back() + 1;

	return 0;
}

int pereads_cluster::project_chain(const map<PI32, PI32> &jm, vector<int32_t> &v)
{
	assert(v.size() % 2 == 0);
	if(v.size() <= 0) return 0;

	for(int k = 0; k < v.size() / 2; k++)
	{
		PI32 p(v[k * 2 + 0], v[k * 2 + 1]);
		map<PI32, PI32>::const_iterator it = jm.find(p);
		if(it == jm.end()) continue;
		v[k * 2 + 0] = it->second.first;
		v[k * 2 + 1] = it->second.second;
	}
	return 0;
}
