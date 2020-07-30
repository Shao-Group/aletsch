/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#include <cassert>
#include <cstdio>
#include <cmath>
#include <climits>

#include "constants.h"
#include "bundle.h"
#include "util.h"

bundle::bundle()
{
	tid = -1;
	chrm = "";
	lpos = 1 << 30;
	rpos = 0;
	strand = '.';
}

int bundle::add_hit_intervals(const hit &ht, bam1_t *b)
{
	add_hit(ht);
	add_intervals(b);
	return 0;
}

int bundle::add_hit(const hit &ht)
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
	if(hits.size() <= 1) strand = ht.strand;
	assert(strand == ht.strand);

	return 0;
}

int bundle::add_intervals(bam1_t *b)
{
	int32_t p = b->core.pos;
	uint32_t *cigar = bam_get_cigar(b);

    for(int k = 0; k < b->core.n_cigar; k++)
	{
		if (bam_cigar_type(bam_cigar_op(cigar[k]))&2)
		{
			p += bam_cigar_oplen(cigar[k]);
		}

		if(bam_cigar_op(cigar[k]) == BAM_CMATCH)
		{
			int32_t s = p - bam_cigar_oplen(cigar[k]);
			mmap += make_pair(ROI(s, p), 1);
		}

		if(bam_cigar_op(cigar[k]) == BAM_CINS)
		{
			imap += make_pair(ROI(p - 1, p + 1), 1);
		}

		if(bam_cigar_op(cigar[k]) == BAM_CDEL)
		{
			int32_t s = p - bam_cigar_oplen(cigar[k]);
			imap += make_pair(ROI(s, p), 1);
		}
	}
	return 0;
}

bool bundle::overlap(const hit &ht) const
{
	if(mmap.find(ROI(ht.pos, ht.pos + 1)) != mmap.end()) return true;
	if(mmap.find(ROI(ht.rpos - 1, ht.rpos)) != mmap.end()) return true;
	return false;
}

int bundle::clear()
{
	tid = -1;
	chrm = "";
	lpos = 1 << 30;
	rpos = 0;
	strand = '.';
	hits.clear();
	mmap.clear();
	imap.clear();
	junctions.clear();
	return 0;
}

int bundle::compute_strand(int libtype)
{
	if(libtype != UNSTRANDED) assert(strand != '.');
	if(libtype != UNSTRANDED) return 0;

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

int bundle::check_left_ascending()
{
	for(int i = 1; i < hits.size(); i++)
	{
		int32_t p1 = hits[i - 1].pos;
		int32_t p2 = hits[i].pos;
		assert(p1 <= p2);
	}
	return 0;
}

int bundle::check_right_ascending()
{
	for(int i = 1; i < hits.size(); i++)
	{
		int32_t p1 = hits[i - 1].rpos;
		int32_t p2 = hits[i].rpos;
		assert(p1 <= p2);
	}
	return 0;
}

int bundle::print(int index)
{
	printf("bundle%d: ", index);

	// statistic xs
	int n0 = 0, np = 0, nq = 0;
	int spliced = 0;
	for(int i = 0; i < hits.size(); i++)
	{
		if(hits[i].xs == '.') n0++;
		if(hits[i].xs == '+') np++;
		if(hits[i].xs == '-') nq++;
		if(hits[i].spos.size() >= 1) spliced++;
	}

	printf("tid = %d, range = %s:%d-%d, orient = %c, #hits = %lu, #spliced = %d, +/-/. = %d / %d / %d\n",
			tid, chrm.c_str(), lpos, rpos, strand, hits.size(), spliced, np, nq, n0);

	return 0;
}

int bundle::build_junctions()
{
	map<PI32, vector<int> > m;
	for(int i = 0; i < bd.hits.size(); i++)
	{
		vector<int32_t> v = bd.hits[i].spos;
		if(v.size() == 0) continue;

		for(int k = 0; k < v.size() / 2; k++)
		{
			PI p(v[k * 2 + 0], v[k * 2 + 1]);
			if(m.find(p) == m.end())
			{
				vector<int> hv;
				hv.push_back(i);
				m.insert(pair< PI32, vector<int> >(p, hv));
			}
			else
			{
				m[p].push_back(i);
			}
		}
	}

	map<PI32, vector<int> >::iterator it;
	for(it = m.begin(); it != m.end(); it++)
	{
		vector<int> &v = it->second;

		int32_t p1 = it->first.first;
		int32_t p2 = it->first.second;

		junction jc(p1, p2, v.size());

		if(jc.count < cfg.min_junction_support) continue;

		for(int k = 0; k < v.size(); k++)
		{
			const hit &h = bd.hits[v[k]];
			jc.nm += h.nm;
			if(h.xs == '.') jc.xs0++;
			if(h.xs == '+') jc.xs1++;
			if(h.xs == '-') jc.xs2++;
		}

		//printf("junction: %s:%d-%d (%d, %d, %d) %d\n", chrm.c_str(), p1, p2, s0, s1, s2, s1 < s2 ? s1 : s2);

		if(jc.xs1 > jc.xs2) jc.strand = '+';
		else if(jc.xs1 < jc.xs2) jc.strand = '-';
		else jc.strand = '.';
		junctions.push_back(jc);

		/*
		uint32_t max_qual = 0;
		for(int k = 0; k < v.size(); k++)
		{
			hit &h = bd.hits[v[k]];
			if(h.qual > max_qual) max_qual = h.qual;
		}
		assert(max_qual >= min_max_boundary_quality);
		*/
	}
	return 0;
}
