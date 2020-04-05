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
	for(MVII::const_iterator x = pmap.begin(); x != pmap.end(); x++)
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
	return 0;
}
