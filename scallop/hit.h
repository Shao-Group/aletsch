/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#ifndef __HIT_H__
#define __HIT_H__

#include <string>
#include <vector>

#include "hit_base.h"

using namespace std;

class hit: public hit_base
{
public:
	//hit(int32_t p);
	hit(bam1_t *b);
	hit(const hit &h);
	~hit();
	hit& operator=(const hit &h);

public:
	vector<int64_t> itvm;					// matched interval
	vector<int64_t> itvi;					// insert interval
	vector<int64_t> itvd;					// delete interval

public:
	int set_intervals(bam1_t *b);
};

#endif
