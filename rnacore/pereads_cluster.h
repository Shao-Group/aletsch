/*
Part of meta-scallop
(c) 2020 by Mingfu Shao, The Pennsylvania State University
See LICENSE for licensing.
*/

#ifndef __PEREADS_CLUSTER_H__
#define __PEREADS_CLUSTER_H__

#include "hit.h"
#include <cstdint>
#include <vector>

using namespace std;

class pereads_cluster
{
public:
	pereads_cluster();

public:
	vector<int32_t> chain1;			// list of intron-chain coordinates
	vector<int32_t> chain2;			// list of intron-chain coordinates
	vector<int32_t> bounds;			// lpos1, rpos1, lpos2, rpos2
	vector<int32_t> extend;			// lexon1, rexon1, lexon2, rexon2
	vector<hit> hits1;				// hits in this cluster (when bridged reads needes to be reported)
	vector<hit> hits2;				// hits in this cluster (when bridged reads needes to be reported)
	int count;						// number of such reads in this cluster

public:
	int clear();
	int print(int k) const;
};

#endif
