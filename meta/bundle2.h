/*
Part of aletsch
(c) 2020 by Mingfu Shao, The Pennsylvania State University
See LICENSE for licensing.
*/

#ifndef __COMBINED_BUNDLE_H__
#define __COMBINED_BUNDLE_H__

#include <vector>
#include <stdint.h>
#include <iostream>
#include <string>

#include "bundle.h"
#include "parameters.h"
#include "sample_profile.h"

using namespace std;

class bundle2 : public bundle
{
public:
	bundle2(const parameters &cfg);
	bundle2(const parameters &cfg, bundle &&bd, int sid);
	bundle2(const bundle2 &cb) = default;
	bundle2(bundle2 &&cb) = default;

public:
	const parameters &cfg;
	int sid;
	string gid;
	int num_combined;

public:
	int bridge(const sample_profile &sp);
	int print(int index);
	int clear();
};

#endif
