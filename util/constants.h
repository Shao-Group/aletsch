/*
Part of Scallop Transcript Assembler
Part of meta-scallop 
(c) 2017 by Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
(c) 2020 by Mingfu Shao, The Pennsylvania State University.
See LICENSE for licensing.
*/

#ifndef __CONSTANTS_H__
#define __CONSTANTS_H__

#include <stdint.h>
#include <map>
#include <set>
#include <sstream>

using namespace std;

#define START_BOUNDARY 1
#define END_BOUNDARY 2
#define LEFT_SPLICE 3
#define RIGHT_SPLICE 4
#define LEFT_RIGHT_SPLICE 5
#define MIDDLE_CUT 6

#define TRIVIAL 0
#define NORMAL 1

// five types for decomposition
#define SMALLEST_EDGE 0
#define NEGLIGIBLE_EDGE 1
#define SPLITTABLE_SIMPLE 2
#define SPLITTABLE_HYPER 3
#define UNSPLITTABLE_SINGLE 4
#define UNSPLITTABLE_MULTIPLE 5
#define TRIVIAL_VERTEX 6

#define EMPTY -1
#define UNSTRANDED 0
#define FR_FIRST 1
#define FR_SECOND 2

// positions
#define IDENTICAL 0
#define FALL_RIGHT 1
#define FALL_LEFT 2
#define CONTAINED 3
#define CONTAINING 4
#define EXTEND_RIGHT 5
#define EXTEND_LEFT 6
#define NESTED 7
#define NESTING 8
#define CONFLICTING 9

// how to combine transcript
#define ADD_TRANSCRIPT_COVERAGE_SUM 1
#define ADD_TRANSCRIPT_COVERAGE_MAX 2

const static string position_names[] = {"identical", "fall-right", "fall-left", "contained", "containing", "extend_right", "extend_left", "nested", "nesting", "conflicting"};

typedef pair<int, int> PI;
typedef pair<double, int> DI;
typedef pair<int32_t, int32_t> PI32;
typedef pair<int32_t, DI> PIDI;
typedef pair<PI32, DI> PPDI;
typedef pair<PI, double> PPID;
typedef pair<int32_t, set<int> > PISI;
typedef map<int, double> MID;
typedef map<int32_t, set<int> > MISI;

#endif
