/*
Part of meta-scallop
(c) 2020 by Mingfu Shao, The Pennsylvania State University
See LICENSE for licensing.
*/

#ifndef __SAMPLE_PROFILE_H__
#define __SAMPLE_PROFILE_H__

#include <vector>
#include <string>

using namespace std;

class sample_profile
{
public:
	int sample_id;
	string file_name;
	int library_type;
	double insertsize_ave;
	double insertsize_std;
	int insertsize_low;
	int insertsize_high;
	int insertsize_median;
	vector<double> insertsize_profile;
};

#endif
