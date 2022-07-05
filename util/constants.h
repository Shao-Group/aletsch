/*
Part of aletsch 
(c) 2020 by Mingfu Shao, The Pennsylvania State University.
Part of Scallop Transcript Assembler
(c) 2017 by Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#ifndef __CONSTANTS_H__
#define __CONSTANTS_H__

#include <stdint.h>
#include <array>
#include <map>
#include <set>
#include <sstream>
#include <vector>

using namespace std;

#define START_BOUNDARY 1
#define END_BOUNDARY 2
#define LEFT_SPLICE 3
#define RIGHT_SPLICE 4
#define LEFT_RIGHT_SPLICE 5
#define MIDDLE_CUT 6
#define LEFT_MIXED 7
#define RIGHT_MIXED 8

#define TRIVIAL 0
#define NORMAL 1

// five types for decomposition
#define SMALLEST_EDGE 0
#define NEGLIGIBLE_EDGE 1
#define SPLITTABLE_SIMPLE 2
#define SPLITTABLE_HYPER 3
#define SPLITTABLE_PURE 4
#define UNSPLITTABLE_SINGLE 5
#define UNSPLITTABLE_MULTIPLE 6
#define TRIVIAL_VERTEX 7
#define MIXED_DIVIDED 8
#define MIXED_BLOCKED 9
#define MIXED_TRIVIAL 10
#define MIXED_TANGLED 11
#define MIXED_SPLITTABLE 12

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
#define TRANSCRIPT_COUNT_ADD_COVERAGE_ADD 1
#define TRANSCRIPT_COUNT_ADD_COVERAGE_NUL 2
#define TRANSCRIPT_COUNT_MAX_COVERAGE_MAX 3

// RNA-seq data type
#define NUM_DATA_TYPES 6
#define DEFAULT 0
#define PAIRED_END 1
#define SINGLE_END 2
#define PACBIO_CCS 3
#define PACBIO_SUB 4
#define ONT 5

const static string position_names[] = {"identical", "fall-right", "fall-left", "contained", "containing", "extend_right", "extend_left", "nested", "nesting", "conflicting"};

typedef pair<int, int> PI;
typedef pair<double, int> DI;
typedef pair<int32_t, int32_t> PI32;
typedef pair<PI32, int> TI32;
typedef pair<int32_t, DI> PIDI;
typedef pair<PI32, DI> PPDI;
typedef pair<TI32, DI> PTDI;
typedef pair<PI, double> PPID;
typedef pair<int32_t, set<int> > PISI;
typedef map<int, double> MID;
typedef map<int32_t, set<int> > MISI;
typedef array<int, 3> AI3;
typedef pair<vector<int32_t>, int> PVI;
typedef pair<vector<int32_t>, AI3> PVI3;

#endif
