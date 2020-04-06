/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <iostream>
#include <ctime>
#include <cassert>
#include <sstream>

#include "parameters.h"
#include "incubator.h"

using namespace std;

int main(int argc, const char **argv)
{
	srand(time(0));	// TODO
	parameters cfg;

	if(argc == 1)
	{
		cfg.print_copyright();
		cfg.print_help();
		printf("\n");
		cfg.print_logo();
		return 0;
	}

	cfg.parse_arguments(argc, argv);
	cfg.print_command_line(argc, argv);

	incubator icbt(cfg);
	icbt.resolve();

	return 0;
}
