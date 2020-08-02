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
		printf("error: adding empty chain with id %d to chain_set\n", h);
		return 0;
	}

	if(h >= 0 && hmap.find(h) != hmap.end()) 
	{
		printf("error: id %d has already been added to chain_set\n", h);
		return 0;
	}

	int32_t p = v[0];
	if(pmap.find(p) == pmap.end())
	{
		vector<PVI> vv;
		vv.push_back(PVI(v, 1));
		chains.push_back(vv);
		int n = chains.size() - 1;
		pmap.insert(make_pair(p, n));
		if(h >= 0) hmap.insert(make_pair(h, PI(n, 0)));
	}
	else
	{
		int k = pmap[p];
		assert(k >= 0 && k < chains.size());
		vector<PVI> &vv = chains[k];
		bool found = false;
		for(int i = 0; i < vv.size(); i++)
		{
			if(vv[i].first == v)
			{
				if(h >= 0) hmap.insert(make_pair(h, PI(k, i)));
				vv[i].second++;
				found = true;
				break;
			}
		}

		if(found == false)
		{
			vv.push_back(PVI(v, 1));
			int n = vv.size() - 1;
			if(h >= 0) hmap.insert(make_pair(h, PI(k, n)));
		}
	}
	return 0;
}

int chain_set::remove(int h)
{
	if(hmap.find(h) == hmap.end()) return 0;
	PI p = hmap[h];
	assert(p.first >= 0 && p.first < chains.size());
	assert(p.second >= 0 && p.second < chains[p.first].size());
	chains[p.first][p.second].second--;
	if(chains[p.first][p.second].second <= 0) chains[p.first][p.second].second = 0;
	hmap.erase(h);
	return 0;
}

PVI chain_set::get(int h) const
{
	PVI pvi;
	pvi.second = -1;
	if(h < 0) return pvi;
	if(hmap.find(h) == hmap.end()) return pvi;
	const auto &p = hmap.at(h);
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

int chain_set::clear()
{
	hmap.clear();
	pmap.clear();
	chains.clear();
	return 0;
}

vector<int32_t> chain_set::get_splices()
{
	set<int32_t> s;
	for(int i = 0; i < chains.size(); i++)
	{
		for(int j = 0; j < chains[i].size(); j++)
		{
			PVI &p = chains[i][j];
			if(p.second <= 0) continue;
			for(int k = 0; k < p.first.size(); k++)
			{
				s.insert(p.first[k]);
			}
		}
	}

	vector<int32_t> v(s.begin(), s.end());
	sort(v.begin(), v.end());
	return v;
}
