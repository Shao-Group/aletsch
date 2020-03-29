#ifndef __HIT_SET_H__
#define __HIT_SET_H__

#include "rcluster.h"
#include "hit.h"

using namespace std;

class hit_set
{
public:
	hit_set(vector<hit> &hits);

public:
	vector<hit> &hits;			

	vector<bool> paired;			
	vector<PRC> vprc;
	vector<rcluster> vsrc;			

public:
	int build_paired_reads_clusters();
};

#endif
