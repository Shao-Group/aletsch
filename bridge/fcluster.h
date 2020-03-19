/*
Part of Coral
(c) 2019 by Mingfu Shao and The Pennsylvania State University.
See LICENSE for licensing.
*/

#ifndef __FRAGMENT_CLUSTER_H__
#define __FRAGMENT_CLUSTER_H__

#include <vector>
#include "phase.h"
#include "fragment.h"

using namespace std;

class fcluster
{
public:
	vector<fragment> frset;			// set of fragments in this cluster
	vector<int> v1;					// vlist for mate1
	vector<int> v2;					// vlist for mate2
	phase bbp;						// best briding path

public:
	int clear();
	int print(int k) const;
	const vector<int> & get_vlist() const;
};

bool compare_fcluster_v1_v2(const fcluster &fx, const fcluster &fy);

#endif
