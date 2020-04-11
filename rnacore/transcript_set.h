#ifndef __TRANSCRIPT_SET_H__
#define __TRANScRIPT_SET_H__

#include <map>
#include <vector>

#include "transcript.h"

using namespace std;

class transcript_set
{
public:
	map<size_t, vector<transcript>> mt;

public:
	int add(transcript &t, int count, int mode = 1);
	int add(const transcript &t, int mode = 1);
	int add(const vector<transcript> &t, int mode = 1);
	int add(const transcript_set &ts, int mode = 1);
	int add_duplicates(const transcript_set &ts, int mode = 1);
	bool query(const transcript &t) const;
	vector<transcript> get_duplicate_transcripts() const;
	vector<transcript> get_transcripts() const;
};

#endif
