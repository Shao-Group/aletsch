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
	parameters cfg;

	if(argc == 1)
	{
		cfg.print_copyright();
		cfg.print_help();
		printf("\n");
		return 0;
	}

	cfg.parse_arguments(argc, argv);
	cfg.print_command_line(argc, argv);

	incubator icbt(cfg);
	icbt.resolve();

	return 0;
}
