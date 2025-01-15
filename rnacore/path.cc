/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#include "path.h"

#include <cassert>
#include <cstdio>

path::path()
{
	v.clear();
    junc.clear();
	abd = 0;
    weight = 0;
    conf = 0;
	reads = 0;
	length = 0;
	strand = '.';
    count = 0;
	id = "";
}

path::~path()
{}

int path::clear()
{
	v.clear();
    junc.clear();
	abd = 0;
    weight = 0;
    conf = 0;
	reads = 0;
	length = 0;
	strand = '.';
    count = 0;
	id = "";
	return 0;
}

int path::print(int index) const
{
	if(v.size() == 0) return 0;
	printf("path %d: weight = %.2lf, abundance = %.2lf, confidence = %.2lf, count = %d, length = %d, ave-reads = %.2lf, strand = %c, vertices = ", index, weight, abd, conf, count, length, reads / length, strand);
	for(int i = 0; i < v.size() - 1; i++)
	{
		printf("%d, ", v[i]);
	}
	printf("%d, junction = ", v[v.size() - 1]);

    for(int i = 0; i < junc.size(); i++)
    {
        printf("(%d,%d), ", junc[i].first, junc[i].second);
    }
    printf("\n");
	return 0;
}

vector<int> path::index(int n) const
{
	vector<int> vv;
	vv.resize(n, -1);
	for(int i = 1; i < v.size(); i++)
	{
		int s = v[i - 1];
		int t = v[i];
		assert(s >= 0 && s < n);
		assert(t >= 0 && t < n);
		vv[s] = t;
	}
	return vv;
}
