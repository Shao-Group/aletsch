/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
Part of Coral package
(c) 2019 by Mingfu Shao, The Pennsylvania State University
See LICENSE for licensing.
*/

#include <cstring>
#include <cassert>
#include <cstdio>
#include <sstream>
#include <cmath>

#include "hit.h"
#include "config.h"
#include "util.h"

hit& hit::operator=(const hit &h)
{
	hit_base::operator=(h);
	hid = h.hid;
	vlist = h.vlist;
	paired = h.paired;
	bridged = h.bridged;
	qhash = h.qhash;
	next = h.next;
	return *this;
}

hit::hit(const hit &h) 
	:hit_base(h)
{
	hid = h.hid;
	vlist = h.vlist;
	paired = h.paired;
	bridged = h.bridged;
	qhash = h.qhash;
	next = h.next;
}

hit::hit(bam1_t *b, int id) 
	:hit_base(b), hid(id)
{}

int hit::get_aligned_intervals(vector<int64_t> &v) const
{
	v.clear();
	int32_t p1 = pos;
	for(int k = 0; k < spos.size(); k++)
	{
		int32_t p2 = high32(spos[k]);
		v.push_back(pack(p1, p2));
		p1 = low32(spos[k]);
	}
	v.push_back(pack(p1, rpos));
	return 0;
}

string hit::get_qname(bam1_t *b)
{
	char buf[1024];
	char *q = bam_get_qname(b);
	int l = strlen(q);
	memcpy(buf, q, l);
	buf[l] = '\0';
	return string(buf);
}

vector<int> encode_vlist(const vector<int> &v)
{
	vector<int> vv;
	if(v.size() <= 0) return vv;

	int p = v[0];
	int k = 1;
	for(int i = 1; i < v.size(); i++)
	{
		if(v[i] == v[i - 1] + 1)
		{
			k++;
		}
		else
		{
			assert(k >= 1);
			vv.push_back(p);
			vv.push_back(k);
			p = v[i];
			k = 1;
		}
	}
	vv.push_back(p);
	vv.push_back(k);

	/*
	printf("encode: (");
	printv(v);
	printf(") -> (");
	printv(vv);
	printf(")\n");
	*/
	return vv;
}

vector<int> decode_vlist(const vector<int> &v)
{
	vector<int> vv;
	assert(v.size() % 2 == 0);
	if(v.size() <= 0) return vv;

	for(int i = 0; i < v.size() / 2; i++)
	{
		int p = v[i * 2 + 0];
		int k = v[i * 2 + 1];
		for(int j = p; j < p + k; j++)
		{
			vv.push_back(j);
		}
	}
	return vv;
}
