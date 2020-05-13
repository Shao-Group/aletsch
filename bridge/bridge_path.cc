/*
Part of meta-scallop
(c) 2020 by Mingfu Shao, The Pennsylvania State University
See LICENSE for licensing.
*/

#include "bridge_path.h"
#include "util.h"

#include <cassert>
#include <cstdio>

bridge_path::bridge_path()
{
	v.clear();
	chain.clear();
	whole.clear();
	type = 0;
	score = 0;
	count = 0;
	strand = 0;
	stack.clear();
}

bridge_path::~bridge_path()
{}

int bridge_path::clear()
{
	type = 0;
	v.clear();
	chain.clear();
	whole.clear();
	score = 0;
	count = 0;
	strand = 0;
	stack.clear();
	return 0;
}

int bridge_path::print(int index) const
{
	if(v.size() == 0) return 0;
	printf("bridge_path %d: score = %.2lf, count = %d, strand = %d, stack = ( ", index, score, count, strand);
	printv(stack);
	printf("), v = ( ");
	printv(v);
	printf("), chain = ( ");
	printv(chain);
	printf("), whole = ( ");
	printv(whole);
	printf(")\n");
	return 0;
}

bool compare_bridge_path_vertices(const bridge_path &p1, const bridge_path &p2)
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

bool compare_bridge_path_score(const bridge_path &p1, const bridge_path &p2)
{
	if(p1.score > p2.score) return true;
	else return false;
}

bool compare_bridge_path_stack(const bridge_path &p1, const bridge_path &p2)
{
	for(int k = 0; k < p1.stack.size() && k < p2.stack.size(); k++)
	{
		if(p1.stack[k] > p2.stack[k]) return true;
		if(p1.stack[k] < p2.stack[k]) return false;
	}
	return p1.stack.size() > p2.stack.size();
}
