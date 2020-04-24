/*
Part of meta-scallop
(c) 2020 by Mingfu Shao, The Pennsylvania State University
See LICENSE for licensing.
*/


#ifndef __TRANSCRIPT_SET_H__
#define __TRANSCRIPT_SET_H__

#include <map>
#include <vector>
#include <string>

#include "transcript.h"

using namespace std;

class trans_item
{
public:
	transcript trst;
	int count;
	set<int> samples;
};

class transcript_set
{
public:
	string chrm;
	char strand;
	map<size_t, vector<trans_item>> mt;

public:
	int add(const trans_item &t, int mode);
	int add(const transcript &t, int count, int sid, int mode);
	int add(const transcript_set &ts, int min_count, int mode);
	int increase_count(int count);
	int print() const;
	vector<transcript> get_transcripts(int min_count) const;
	pair<bool, trans_item> query(const transcript &t) const;
};

#endif
