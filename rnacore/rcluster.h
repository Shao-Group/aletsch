#ifndef __READ_CLUSTER_H__
#define __READ_CLUSTER_H__

#include <cstdint>
#include <vector>

using namespace std;

class rcluster
{
public:
	vector<int32_t> vv;			// list of vertices / intron-chain coordinates
	vector<int32_t> vl;			// list of left offset
	vector<int32_t> vr;			// list of right offset 
	vector<int> cc;				// count of each offset

public:
	int clear();
	int print(int k) const;
};

typedef pair<rcluster, rcluster> PRC;

#endif
