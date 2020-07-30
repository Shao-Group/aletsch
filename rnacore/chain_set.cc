/*
Part of aletsch
(c) 2020 by Mingfu Shao, The Pennsylvania State University
See LICENSE for licensing.
*/

#include "essential.h"
#include "constants.h"
#include "chain_set.h"

int chain_set::add(const vector<int32_t> &v, int h)
{
	if(v.size() <= 0)
	{
		printf("error: adding empty chain with id to chain_set\n", h);
		return 0;
	}

	if(hmap.find(h) != h.end()) 
	{
		printf("error: id %d has already been added to chain_set\n", h);
		return 0;
	}

	int32_t p = v[0];
	if(pmap.find(p) == pmap.end())
	{
		vector<vector<int32_t>> vv;
		vv.push_back(v);
		chains.push_back(vv);
		int n = chains.size() - 1;
		pmap.insert(make_pair(p, n));
		hmap.insert(make_pair(h, PI(n, 0)));
	}
	else
	{
		int k = pmap[p];
		assert(k >= 0 && k < chains.size());
		vector<vector<int32_t>> &vv = chains[k];
		bool found = false;
		for(int i = 0; i < vv.size(); i++)
		{
			if(vv[i] == v)
			{
				hmap.insert(make_pair(h, PI(k, i)));
				found = true;
				break;
			}
		}

		if(found == false)
		{
			vv.push_back(v);
			int n = vv.size() - 1;
			hmap.insert(make_pair(h, PI(k, n)));
		}
	}
	return 0;
}

vector<int32_t> chain_set::get_chain(int h)
{
	vector<int32_t> v;
	if(hmap.find(h) == hmap.end()) return v;

	PI p = hmap[h];
	assert(p.first >= 0 && p.first < chains.size());
	assert(p.second >= 0 && p.second < chains[p.first].size());
	return chains[p.first][p.second];
}

int chain_set::print()
{
	map<int, int> m;
	for(int i = 0; i < chains.size(); i++)
	{
		int n = chains[i].size();
		if(m.find(n) == m.end()) m.insert(make_pair(n, 1));
		else m[n]++;
	}

	printf("chain_set: %lu groups, %lu stored hits\n", pmap.size(), hmap.size());
	for(auto &x : m)
	{
		printf("chain_set: %d groups with %d chains\n", x.second, x.first);
	}
	return 0;
}
