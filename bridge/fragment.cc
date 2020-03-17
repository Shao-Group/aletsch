/*
Part of Coral, an efficient tool to bridge mate pairs
(c) 2018 by Mingfu Shao and The Pennsylvania State University.
See LICENSE for licensing.
*/

#include "fragment.h"
#include "util.h"
#include <cstdio>

fragment::fragment(int x1, int x2)
	: h1(x1), h2(x2)
{
	b1 = false;
	b2 = false;
	k1l = 0;
	k1r = 0;
	k2l = 0;
	k2r = 0;
	lpos = -1;
	rpos = -1;
	cnt = 1;
}

int fragment::print(int index)
{
	printf("fragment %d: (%d, %d), cnt = %d, lpos = %d, rpos = %d, len = %d, b1 = %c, b2 = %c\n",
			index, h1, h2, cnt, lpos, rpos, rpos - lpos, b1 ? 'T' : 'F', b2 ? 'T' : 'F');
	return 0;
}

int fragment::clear()
{
	h1 = -1;
	h2 = -1;
	b1 = false;
	b2 = false;
	k1l = 0;
	k1r = 0;
	k2l = 0;
	k2r = 0;
	lpos = -1;
	rpos = -1;
	cnt = 1;
	return 0;
}

bool fragment::equal(const fragment &f) const
{
	/*
	if(h1->vlist != f.h1->vlist) return false;
	if(h2->vlist != f.h2->vlist) return false;
	*/
	if(lpos != f.lpos) return false;
	if(rpos != f.rpos) return false;
	if(k1l != f.k1l) return false;
	if(k1r != f.k1r) return false;
	if(k2l != f.k2l) return false;
	if(k2r != f.k2r) return false;
	if(b1 != f.b1) return false;
	if(b2 != f.b2) return false;
	return true;
}
