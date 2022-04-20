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
#include "bundle_base.h"
#include "sample_profile.h"

using namespace std;

class bundle : public bundle_base
{
public:
	bundle(const parameters &cfg, const sample_profile &sp);
	bundle(const parameters &cfg, const sample_profile &sp, bundle_base &&bb);
	bundle(const bundle &cb) = default;
	bundle(bundle &&cb) = default;

public:
	const parameters &cfg;
	const sample_profile &sp;
	int num_combined;
	string gid;

public:
	int set_gid(int instance, int subindex);
	int copy_meta_information(const bundle &bb);
	int combine(const bundle &bb);
	int bridge();
	int clear();
};

#endif
