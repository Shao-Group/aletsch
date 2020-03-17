/*
Part of Coral, an efficient tool to bridge mate pairs
(c) 2018 by Mingfu Shao and The Pennsylvania State University.
See LICENSE for licensing.
*/

#ifndef __PIER_H__
#define __PIER_H__

#include <vector>
#include <stdint.h>

#include "phase.h"

using namespace std;

class pier
{
public:
	pier(int s, int t);

public:
	bool operator<(const pier &p) const;

public:
	int bs;					// source of bridge
	int bt;					// target of bridge
	vector<phase> phases;	// bridging candidates 
};

#endif
