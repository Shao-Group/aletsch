/*
Part of meta-scallop
(c) 2020 by Mingfu Shao, The Pennsylvania State University
See LICENSE for licensing.
*/

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
