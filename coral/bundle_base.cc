/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#include <cassert>
#include <cstdio>
#include <cmath>
#include <climits>

#include "bundle_base.h"
#include "config.h"

bundle_base::bundle_base()
{
	tid = -1;
	chrm = "";
	lpos = 1 << 30;
	rpos = 0;
	strand = '.';
}

bundle_base::~bundle_base()
{}

int bundle_base::add_hit(const hit &ht)
{
	// store new hit
	hits.push_back(ht);

	// calcuate the boundaries on reference
	if(ht.pos < lpos) lpos = ht.pos;
	if(ht.rpos > rpos) rpos = ht.rpos;

	// set tid
	if(tid == -1) tid = ht.tid;
	assert(tid == ht.tid);

	// set strand
	//if(hits.size() <= 1) strand = ht.strand;
	assert(strand == ht.strand || ht.strand == '.');

	// DEBUG
	/*
	if(strand != ht.strand)
	{
		printf("strand = %c, ht.strand = %c, ht.xs = %c,\n", strand, ht.strand, ht.xs);
	}
	*/

	// TODO: do not add to interval-map (to reduce memory; for coral)
	/*
	for(int k = 0; k < ht.itvm.size(); k++)
	{
		int32_t s = high32(ht.itvm[k]);
		int32_t t = low32(ht.itvm[k]);
		//printf(" add interval %d-%d\n", s, t);
		mmap += make_pair(ROI(s, t), 1);
	}

	for(int k = 0; k < ht.itvi.size(); k++)
	{
		int32_t s = high32(ht.itvi[k]);
		int32_t t = low32(ht.itvi[k]);
		imap += make_pair(ROI(s, t), 1);
	}

	for(int k = 0; k < ht.itvd.size(); k++)
	{
		int32_t s = high32(ht.itvd[k]);
		int32_t t = low32(ht.itvd[k]);
		imap += make_pair(ROI(s, t), 1);
	}
	*/

	return 0;
}

bool bundle_base::overlap(const hit &ht) const
{
	if(mmap.find(ROI(ht.pos, ht.pos + 1)) != mmap.end()) return true;
	if(mmap.find(ROI(ht.rpos - 1, ht.rpos)) != mmap.end()) return true;
	return false;
}

int bundle_base::clear()
{
	tid = -1;
	chrm = "";
	lpos = 1 << 30;
	rpos = 0;
	strand = '.';
	hits.clear();
	mmap.clear();
	imap.clear();
	return 0;
}

int bundle_base::compute_strand()
{
	if(library_type != UNSTRANDED) assert(strand != '.');
	if(library_type != UNSTRANDED) return 0;

	int n0 = 0, np = 0, nq = 0;
	for(int i = 0; i < hits.size(); i++)
	{
		if(hits[i].xs == '.') n0++;
		if(hits[i].xs == '+') np++;
		if(hits[i].xs == '-') nq++;
	}

	if(np > nq) strand = '+';
	else if(np < nq) strand = '-';
	else strand = '.';

	return 0;
}

int bundle_base::check_left_ascending()
{
	for(int i = 1; i < hits.size(); i++)
	{
		int32_t p1 = hits[i - 1].pos;
		int32_t p2 = hits[i].pos;
		assert(p1 <= p2);
	}
	return 0;
}

int bundle_base::check_right_ascending()
{
	for(int i = 1; i < hits.size(); i++)
	{
		int32_t p1 = hits[i - 1].rpos;
		int32_t p2 = hits[i].rpos;
		assert(p1 <= p2);
	}
	return 0;
}
