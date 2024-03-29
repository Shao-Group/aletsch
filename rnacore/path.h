/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#ifndef __PATH_H__
#define __PATH_H__

#include <vector>

using namespace std;

class path
{
public:
	path();
	~path();

public:
	vector<int> v;
    vector<pair<int, int>> junc;
	int length;
	double abd;
    double weight;
    double conf;
	double reads;
	char strand;
    int count;

public:
	int clear();
	int print(int index) const;
	vector<int> index(int n) const;
};

#endif
