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

	// add splices
	vector<int32_t> v = extract_splices(b);
	hcst.add(v, hits.size() - 1);

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
	hcst.clear();
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

int bundle::build_fragments()
{
	frgs.clear();
	brdg.clear();
	if(hits.size() == 0) return 0;

	int max_index = hits.size() + 1;
	if(max_index > 1000000) max_index = 1000000;

	vector<bool> paired(hits.size(), false);
	vector< vector<int> > vv;
	vv.resize(max_index);

	// first build index
	for(int i = 0; i < hits.size(); i++)
	{
		const hit &h = hits[i];
		// do not use hi; as long as qname, pos and isize are identical
		int k = (h.get_qhash() % max_index + h.pos % max_index + (0 - h.isize) % max_index) % max_index;
		vv[k].push_back(i);
	}

	for(int i = 0; i < hits.size(); i++)
	{
		const hit &h = hits[i];
		if(paired[i] == true) continue;

		int k = (h.get_qhash() % max_index + h.mpos % max_index + h.isize % max_index) % max_index;
		int x = -1;
		for(int j = 0; j < vv[k].size(); j++)
		{
			int u = vv[k][j];
			const hit &z = hits[u];
			if(u == i) continue;
			//if(z.hi != h.hi) continue;
			if(paired[u] == true) continue;
			if(z.pos != h.mpos) continue;
			if(z.isize + h.isize != 0) continue;
			//if(z.qhash != h.qhash) continue;
			if(z.qname != h.qname) continue;
			x = u;
			break;
		}

		if(x == -1) continue;

		assert(i != x);
		//frgs.push_back(PI(hits[i].hid, hits[x].hid));
		frgs.push_back(PI(i, x));
		brdg.push_back(false);

		paired[i] = true;
		paired[x] = true;
	}

	//printf("total hits = %lu, total fragments = %lu\n", hits.size(), frgs.size());
	return 0;
}

int bundle::build_junctions()
{
	map<PI32, vector<int> > m;
	for(int i = 0; i < bd.hits.size(); i++)
	{
		//vector<int32_t> v = bd.hits[i].spos;
		vector<int32_t> v = hcst.get_chain(hits[i].hid);
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
		jcns.push_back(jc);

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

int bundle::update_bridged_fragments(const vector<int> &frlist, const vector<int32_t> &chain, const vector<int32_t> &whole)
{
	return 0;
}
