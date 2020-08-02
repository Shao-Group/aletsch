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

#include "parameters.h"

using namespace std;

class combined_bundle : public bundle
{
public:
	combined_bundle(const parameters &cfg);
	combined_bundle(const combined_bundle &cb) = default;
	combined_bundle(combined_bundle &&cb) = default;

public:
	const parameters &cfg;
	int sid;
	string gid;
	int num_combined;

public:
	int print(int index);
	int clear();
};

#endif
