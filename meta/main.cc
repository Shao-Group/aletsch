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

#include "scallop/config.h"
#include "meta_config.h"
#include "incubator.h"

using namespace std;

int main(int argc, const char **argv)
{
	srand(time(0));

	if(argc == 1)
	{
		print_meta_copyright();
		print_meta_help();
		printf("\n");
		print_meta_logo();
		return 0;
	}

	config cfg;
	cfg.parse_arguments(argc, argv);
	parse_meta_arguments(argc, argv);

	incubator icbt(cfg);
	icbt.resolve();

	return 0;
}
