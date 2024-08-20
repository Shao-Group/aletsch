/*
Part of aletsch
(c) 2024 by Qian Shi, Qimin Zhang, Mingfu Shao, The Pennsylvania State University
See LICENSE for licensing.
*/

#include <cstdio>
#include <iostream>

#include "parameters.h"
#include "incubator.h"
#include "interval_map.h"

using namespace std;

int main(int argc, const char **argv)
{
	//test_join_interval_map();
	//return 0;
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
