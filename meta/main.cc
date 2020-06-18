/*
Part of meta-scallop
(c) 2020 by Mingfu Shao, The Pennsylvania State University
See LICENSE for licensing.
*/

#include <cstdio>
#include <iostream>

#include "parameters.h"
#include "incubator.h"

#include "phase_set.h"

using namespace std;

int test(const phase_set &ps)
{
	MVII::const_iterator it = ps.pmap.begin();
	assert(it == ps.pmap.end());

	for(MVII::const_iterator it = ps.pmap.begin(); it != ps.pmap.end(); it++)
	{
		const vector<int32_t> &v = it->first;
		int c = it->second;
		printv(v);
		printf(", c = %d\n", c);
	}
	return 0;
}

int testps()
{
	phase_set ps;
	vector<int> v = {1, 2, 3, 4};
	ps.add(v, 10);
	ps.add(v, 20);
	//test(ps);
	phase_set ps2;
	ps2.pmap = std::move(ps.pmap);
	//test(ps2);
	test(ps);
	return 0;
}

int main(int argc, const char **argv)
{
	setbuf(stdout, NULL);

	vector<parameters> params(NUM_DATA_TYPES);
	for(int i = 0; i < NUM_DATA_TYPES; i++)
	{
		params[i].set_default(i);
		params[i].parse_arguments(argc, argv, i);
	}

	if(argc == 1)
	{
		params[0].print_copyright();
		params[0].print_help();
		printf("\n");
		return 0;
	}

	params[0].print_command_line(argc, argv);
	incubator icbt(params);
	icbt.resolve();

	return 0;
}
