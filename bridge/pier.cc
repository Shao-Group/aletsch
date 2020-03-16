#include "pier.h"

pier::pier(int s, int t)
	:bs(s), bt(t)
{}

bool pier::operator<(const pier &p) const
{
	if(bs < p.bs) return true;
	if(bs > p.bs) return false;
	return (bt < p.bt);
}
