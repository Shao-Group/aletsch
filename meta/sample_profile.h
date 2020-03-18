#ifndef __SAMPLE_PROFILE_H__
#define __SAMPLE_PROFILE_H__

#include <vector>
using namespace std;

class sample_profile
{
public:
	int library_type;
	double insertsize_ave;
	double insertsize_std;
	int insertsize_low;
	int insertsize_high;
	int insertsize_median;
	vector<double> insertsize_profile;
};

#endif
