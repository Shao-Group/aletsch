/*
Part of meta-scallop
(c) 2020 by Mingfu Shao, The Pennsylvania State University
See LICENSE for licensing.
*/


#include "essential.h"
#include "constants.h"
#include "phase_set.h"

int phase_set::clear()
{
	pmap.clear();
	return 0;
}

int phase_set::add(const vector<int32_t> &v, int c)
{
	if(v.size() <= 0)
	{
		printf("error: adding empty vector to phase_set\n");
		return 0;
	}

	assert(v.size() % 2 == 0);
	if(pmap.find(v) == pmap.end()) pmap.insert(PVII(v, c));
	else pmap[v] += c;
	return 0;
}

int phase_set::print()
{
	for(MVII::iterator it = pmap.begin(); it != pmap.end(); it++)
	{
		const vector<int32_t> &v = it->first;
		int c = it->second;
		printf("phase: count = %d, exons = ( ", c);
		printv(v);
		printf("\n");
	}
	return 0;
}

int phase_set::combine(const phase_set &ps)
{
	for(MVII::const_iterator it = ps.pmap.begin(); it != ps.pmap.end(); it++)
	{
		const vector<int32_t> &v = it->first;
		int c = it->second;
		add(v, c);
	}
	return 0;
}

int phase_set::project_boundaries(const map<int32_t, int32_t> &smap, const map<int32_t, int32_t> &tmap)
{
	phase_set ps;
	for(auto x = pmap.begin(); x != pmap.end(); x++)
	{
		vector<int32_t> v = x->first;
		int c = x->second;
		assert(v.size() % 2 == 0);
		assert(v.size() >= 2);
		map<int32_t, int32_t>::const_iterator is = smap.find(v.front());
		map<int32_t, int32_t>::const_iterator it = tmap.find(v.back());
		if(is != smap.end()) v[0] = is->second;
		if(it != tmap.end()) v[v.size() - 1] = it->second;
		ps.add(v, c);
	}
	pmap = ps.pmap;
	//pmap = std::move(ps.pmap);
	return 0;
}

int phase_set::project_junctions(const map<PI32, PI32> &jm)
{
	phase_set ps;
	for(MVII::iterator x = pmap.begin(); x != pmap.end(); x++)
	{
		const vector<int32_t> &v = x->first;
		int c = x->second;
		assert(v.size() % 2 == 0);
		assert(v.size() >= 2);

		vector<int32_t> vv;
		vv.push_back(v.front());
		for(int k = 0; k < v.size() / 2 - 1; k++)
		{
			PI32 p(v[k * 2 + 1], v[k * 2 + 2]);
			map<PI32, PI32>::const_iterator it = jm.find(p);
			if(it == jm.end())
			{
				vv.push_back(v[k * 2 + 1]);
				vv.push_back(v[k * 2 + 2]);
			}
			else
			{
				vv.push_back(it->second.first);
				vv.push_back(it->second.second);
			}
		}
		vv.push_back(v.back());
		assert(v.size() == vv.size());

		bool succeed = check_increasing_sequence(vv);
		if(succeed == true) ps.add(vv, c);
		//else ps.add(v, c);
	}
	//pmap = std::move(ps.pmap);
	pmap = ps.pmap;
	return 0;
}
