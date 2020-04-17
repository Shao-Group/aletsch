/*
Part of meta-scallop
(c) 2020 by Mingfu Shao, The Pennsylvania State University
See LICENSE for licensing.
*/

#ifndef __PIER_H__
#define __PIER_H__

#include <vector>
#include <stdint.h>

#include "bridge_path.h"

using namespace std;

class pier
{
public:
	pier(int s, int t);

public:
	bool operator<(const pier &p) const;

public:
	int bs;							// source of bridge
	int bt;							// target of bridge
	vector<bridge_path> bridges;	// bridging candidates 
};

#endif
