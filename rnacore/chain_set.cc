/*
Part of aletsch
(c) 2020 by Mingfu Shao, The Pennsylvania State University
See LICENSE for licensing.
*/

#include "essential.h"
#include "constants.h"
#include "chain_set.h"

int chain_set::add(const chain_set &cst)
{
	for(int i = 0; i < cst.chains.size(); i++)
	{
		for(int j = 0; j < cst.chains[i].size(); j++)  
		{
			const PVI3 &p = cst.chains[i][j];
			add(p.first, p.second);
		}
	}
	return 0;
}

int chain_set::add(const vector<int32_t> &v, const AI3 &a)
{
	if(v.size() <= 0)
	{
		printf("error: adding empty chain to chain_set\n");
		return 0;
	}

	int32_t p = v[0];
	if(pmap.find(p) == pmap.end())
	{
		vector<PVI3> vv;
		vv.push_back(PVI3(v, a));
		chains.push_back(vv);
		int n = chains.size() - 1;
		pmap.insert(PI(p, n));
	}
	else
	{
		int k = pmap[p];
		assert(k >= 0 && k < chains.size());
		vector<PVI3> &vv = chains[k];
		bool found = false;
		for(int i = 0; i < vv.size(); i++)
		{
			if(vv[i].first == v)
			{
				vv[i].second[0] += a[0];
				vv[i].second[1] += a[1];
				vv[i].second[2] += a[2];
				found = true;
				break;
			}
		}

		if(found == false) vv.push_back(PVI3(v, a));
	}
	return 0;
}

int chain_set::add(const vector<int32_t> &v, int h, char c)
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

	int xs = 0;
	if(c == '+') xs = 1;
	if(c == '-') xs = 2;

	int32_t p = v[0];
	if(pmap.find(p) == pmap.end())
	{
		AI3 a = {0, 0, 0};
		a[xs] = 1;
		vector<PVI3> vv;
		vv.push_back(PVI3(v, a));
		chains.push_back(vv);
		int n = chains.size() - 1;
		pmap.insert(PI(p, n));
		if(h >= 0) hmap.insert(make_pair(h, AI3({n, 0, xs})));
	}
	else
	{
		int k = pmap[p];
		assert(k >= 0 && k < chains.size());
		vector<PVI3> &vv = chains[k];
		bool found = false;
		for(int i = 0; i < vv.size(); i++)
		{
			if(vv[i].first == v)
			{
				if(h >= 0) hmap.insert(make_pair(h, AI3({k, i, xs})));
				vv[i].second[xs]++;
				found = true;
				break;
			}
		}

		if(found == false)
		{
			AI3 a = {0, 0, 0};
			a[xs] = 1;
			vv.push_back(PVI3(v, a));
			int n = vv.size() - 1;
			if(h >= 0) hmap.insert(make_pair(h, AI3({k, n, xs})));
		}
	}
	return 0;
}

int chain_set::remove(int h)
{
	if(hmap.find(h) == hmap.end()) return 0;
	AI3 p = hmap[h];
	assert(p[0] >= 0 && p[0] < chains.size());
	assert(p[1] >= 0 && p[1] < chains[p[0]].size());
	assert(p[2] >= 0 && p[2] <= 2);
	chains[p[0]][p[1]].second[p[2]]--;
	if(chains[p[0]][p[1]].second[p[2]] <= 0) chains[p[0]][p[1]].second[p[2]] = 0;
	hmap.erase(h);
	return 0;
}

vector<int32_t> chain_set::get_chain(int h) const
{
	vector<int32_t> v;
	if(h < 0) return v;
	if(hmap.find(h) == hmap.end()) return v;
	const auto &p = hmap.at(h);
	return chains[p[0]][p[1]].first;
}

PVI3 chain_set::get(int h) const
{
	PVI3 pvi;
	pvi.second = {-1, -1, -1};
	if(h < 0) return pvi;
	if(hmap.find(h) == hmap.end()) return pvi;
	const auto &p = hmap.at(h);
	return chains[p[0]][p[1]];
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
	map<int, AI3>().swap(hmap);
	map<int32_t, int>().swap(pmap);
	for(int k = 0; k < chains.size(); k++) vector<PVI3>().swap(chains[k]);
	vector<vector<PVI3>>().swap(chains);
	return 0;
}

vector<int32_t> chain_set::get_splices() const
{
	set<int32_t> s;
	for(int i = 0; i < chains.size(); i++)
	{
		for(int j = 0; j < chains[i].size(); j++)
		{
			const PVI3 &p = chains[i][j];
			if(p.second[0] + p.second[1] + p.second[2] <= 0) continue;
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
