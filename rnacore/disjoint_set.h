/*
Part of aletsch
(c) 2020 by Mingfu Shao, The Pennsylvania State University
See LICENSE for licensing.
*/

#ifndef __DISJOINT_SET_H__
#define __DISJOINT_SET_H__

#include "boost/pending/disjoint_sets.hpp"

typedef boost::disjoint_sets<int*, int*> disjoint_set_t;

class disjoint_set
{
public:
	disjoint_set();
	disjoint_set(int n);
	~disjoint_set();

public:
	int init(int n);
	int clear();
	int find_set(int x) {  return ds->find_set(x); }
	int get_size(int p) { return sizes[p]; }
	void set_size(int p, int m) { sizes[p] = m; }
	void link(int x, int y) { ds->link(x, y); }

private:
	int num;
	int *sizes;
	int *ranks;
	int *parents;
	disjoint_set_t *ds;
};

#endif
