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
