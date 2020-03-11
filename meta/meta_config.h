/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#ifndef __META_CONFIG_H__
#define __META_CONFIG_H__

#include "util.h"
#include <stdint.h>
#include <map>
#include <sstream>

using namespace std;

extern int min_supporting_samples;
extern int min_splicing_count;
extern int min_phasing_count;
extern int32_t max_group_boundary_distance;
extern double max_group_junction_distance;
extern bool merge_intersection;
extern int max_threads;
extern int max_combined;

extern int meta_verbose;
extern string meta_version;

// parse arguments
int parse_meta_arguments(int argc, const char ** argv);
int print_meta_copyright();
int print_meta_logo();
int print_meta_help();

#endif
