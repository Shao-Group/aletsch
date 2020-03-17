/*
Part of Coral, an efficient tool to bridge mate pairs
(c) 2018 by Mingfu Shao and The Pennsylvania State University.
See LICENSE for licensing.
*/

#ifndef __FRAGMENT_H__
#define __FRAGMENT_H__

#include <vector>
#include <stdint.h>

using namespace std;

class fragment
{
public:
	fragment(int h1, int h2);

public:
	int h1;			// first mate
	int h2;			// second mate
	int cnt;			// count of the equal hits
	int32_t lpos;		// equals to hits[k1].pos
	int32_t rpos;		// equals to hits[k2].rpos
	int32_t k1l;		// k1-left-outside length
	int32_t k1r;		// k1-right-outside length
	int32_t k2l;		// k2-left-outside length
	int32_t k2r;		// k2-right-outside length
	bool b1;			// whether left mate can be shorten
	bool b2;			// whether right mate can be shorten

public:
	bool equal(const fragment &f) const;
	int print(int index);
	int clear();
};

#endif
