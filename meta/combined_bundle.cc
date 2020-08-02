/*
Part of aletsch
(c) 2020 by Mingfu Shao, The Pennsylvania State University
See LICENSE for licensing.
*/

#include "bridge_solver.h"
#include "parameters.h"
#include "combined_bundle.h"
#include "config.h"
#include "essential.h"
#include <sstream>
#include <algorithm>

combined_bundle::combined_bundle(const parameters &c)
	: cfg(c)
{
	sid = -1;
	num_combined = 0;
}
