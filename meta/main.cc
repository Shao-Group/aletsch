/*
Part of meta-scallop
(c) 2020 by Mingfu Shao, The Pennsylvania State University
See LICENSE for licensing.
*/

#include <cstdio>
#include <iostream>

#include "parameters.h"
#include "incubator.h"

using namespace std;

int main(int argc, const char **argv)
{
	parameters params[NUM_DATA_TYPES];
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
	incubator icbt(params[0]);
	icbt.resolve();

	return 0;
}
