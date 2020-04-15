#ifndef __TRANSCRIPT_SET_H__
#define __TRANSCRIPT_SET_H__

#include <map>
#include <vector>
#include <string>

#include "transcript.h"

using namespace std;

class transcript_set
{
public:
	string chrm;
	char strand;
	map<size_t, vector<transcript>> mt;

public:
	int add(const transcript &t, int mode);
	int add(const transcript_set &ts, int mode);
	int add(const transcript_set &ts, int min_count, int mode);
	int increase_count(int count);
	vector<transcript> get_transcripts(int min_count) const;
	bool query(const transcript &t) const;
	int print() const;
};

#endif
