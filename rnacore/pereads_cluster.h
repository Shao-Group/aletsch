#ifndef __PEREADS_CLUSTER_H__
#define __PEREADS_CLUSTER_H__

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
	int count;						// number of such reads in this cluster

public:
	int clear();
	int print(int k) const;
};

#endif
