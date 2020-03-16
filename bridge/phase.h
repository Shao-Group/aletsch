/*
Part of Scallop Transcript Assembler
(c) 2017 by Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.

Part of Coral, an efficient tool to bridge mate pairs
(c) 2018 by Mingfu Shao and The Pennsylvania State University.

See LICENSE for licensing.
*/

#ifndef __PHASE_H__
#define __PHASE_H__

#include <vector>
#include <stdint.h>

using namespace std;

class phase
{
public:
	phase();
	~phase();

public:
	int type;
	vector<int> v;
	int32_t length;
	double score;
	vector<int> stack;

public:
	int clear();
	int print(int index) const;
	int print_bridge(int index) const;
	vector<int> index(int n) const;
};

bool compare_phase_vertices(const phase &p1, const phase &p2);
bool compare_phase_score(const phase &p1, const phase &p2);
bool compare_phase_stack(const phase &p1, const phase &p2);

#endif
