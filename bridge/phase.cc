/*
Part of Scallop Transcript Assembler
(c) 2017 by Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.

Part of Coral, an efficient tool to bridge mate pairs
(c) 2018 by Mingfu Shao and The Pennsylvania State University.

See LICENSE for licensing.
*/

#include "phase.h"
#include "util.h"

#include <cassert>
#include <cstdio>

phase::phase()
{
	v.clear();
	type = 0;
	score = 0;
	length = 0;
	stack.clear();
}

phase::~phase()
{}

int phase::clear()
{
	type = 0;
	v.clear();
	score = 0;
	length = 0;
	stack.clear();
	return 0;
}

int phase::print(int index) const
{
	if(v.size() == 0) return 0;
	printf("phase %d: score = %.2lf, length = %d, stack = ( ", index, score, length);
	printv(stack);
	printf("), v = ( ");
	printv(v);
	printf("\n");
	return 0;
}

bool compare_phase_vertices(const phase &p1, const phase &p2)
{
	for(int k = 0; k < p1.v.size() && k < p2.v.size(); k++)
	{
		if(p1.v[k] < p2.v[k]) return true;
		if(p1.v[k] > p2.v[k]) return false;
	}
	if(p1.v.size() < p2.v.size()) return true;
	if(p1.v.size() > p2.v.size()) return false;
	return false;
}

bool compare_phase_score(const phase &p1, const phase &p2)
{
	if(p1.score > p2.score) return true;
	else return false;
}

bool compare_phase_stack(const phase &p1, const phase &p2)
{
	for(int k = 0; k < p1.stack.size() && k < p2.stack.size(); k++)
	{
		if(p1.stack[k] > p2.stack[k]) return true;
		if(p1.stack[k] < p2.stack[k]) return false;
	}
	return p1.stack.size() > p2.stack.size();
}
