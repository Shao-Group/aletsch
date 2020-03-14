/*
Part of Scallop Transcript Assembler
(c) 2017 by Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.

Part of Coral, an efficient tool to bridge mate pairs
(c) 2018 by Mingfu Shao and The Pennsylvania State University.

See LICENSE for licensing.
*/

#ifndef __PATH_H__
#define __PATH_H__

#include <vector>
#include <stdint.h>

#include "region.h"

using namespace std;

class path
{
public:
	path();
	~path();

public:
	bool operator< (const path& p) const;

public:
	int type;
	int ex1;
	int ex2;
	vector<int> v;
	vector<int32_t> acc;
	int32_t length;
	int fcindex;
	double abd;
	double prlen;
	double score;
	double reads;

public:
	int clear();
	int print(int index) const;
	int print_bridge(int index) const;
	vector<int> index(int n) const;
	int build_accumulate_length(const vector<region> &regions);
};

bool compare_path_abundance(const path &p1, const path &p2);
bool compare_path_vertices(const path &p1, const path &p2);
bool compare_path_score(const path &p1, const path &p2);

#endif
