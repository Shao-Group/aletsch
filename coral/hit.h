/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
Part of Coral package
(c) 2019 by Mingfu Shao, The Pennsylvania State University
See LICENSE for licensing.
*/

#ifndef __HIT_H__
#define __HIT_H__

#include <string>
#include <vector>

#include "hit_base.h"

using namespace std;

class hit : public hit_base
{
public:
	//hit(int32_t p);
	hit(bam1_t *b, int id);
	hit(const hit &h);
	hit& operator=(const hit &h);

public:
	int hid;								// hit-id
	vector<int> vlist;						// list of spanned vertices in the junction graph
	size_t qhash;							// hash code for qname
	bool paired;							// whether this hit has been paired
	bool bridged;							// whether this hit has been bridged 
	hit *next;								// next hit that is equivalent with current one

public:
 	static string get_qname(bam1_t *b);
	int get_aligned_intervals(vector<int64_t> &v) const;
};

vector<int> encode_vlist(const vector<int> &v);
vector<int> decode_vlist(const vector<int> &v);

#endif
