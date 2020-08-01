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
#include "essential.h"

bundle::bundle()
{
	tid = -1;
	chrm = "";
	lpos = 1 << 30;
	rpos = 0;
	strand = '.';
}

int bundle::build()
{
	build_fragments();
	filter_secondary_hits();
	build_fragments();
	return 0;
}

int bundle::add_hit_intervals(const hit &ht, bam1_t *b)
{
	add_hit(ht);
	add_intervals(b);
	vector<int32_t> v = ht.extract_splices(b);
	if(v.size() >= 1) hcst.add(v, hits.size() - 1);
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
		if(bam_cigar_type(bam_cigar_op(cigar[k]))&2)
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
	fcst.clear();
	mmap.clear();
	imap.clear();
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
	for(int i = 0; i < hits.size(); i++)
	{
		if(hits[i].xs == '.') n0++;
		if(hits[i].xs == '+') np++;
		if(hits[i].xs == '-') nq++;
	}

	printf("tid = %d, range = %s:%d-%d, orient = %c, #hits = %lu, +/-/. = %d / %d / %d\n",
			tid, chrm.c_str(), lpos, rpos, strand, hits.size(), np, nq, n0);

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
		brdg.push_back(0);

		paired[i] = true;
		paired[x] = true;
	}

	//printf("total hits = %lu, total fragments = %lu\n", hits.size(), frgs.size());
	return 0;
}

int bundle::filter_secondary_hits()
{
	set<string> primary;
	for(int i = 0; i < frgs.size(); i++)
	{
		int h1 = frgs[i].first;
		int h2 = frgs[i].second;
		assert(hits[h1].qname == hits[h2].qname);
		if((hits[h1].flag & 0x100) <= 0 && (hits[h2].flag & 0x100) <= 0)
		{
			primary.insert(hits[h1].qname);
		}
	}

	int cnt = 0;
	vector<bool> redundant(hits.size(), false);
	for(int i = 0; i < frgs.size(); i++)
	{
		int h1 = frgs[i].first;
		int h2 = frgs[i].second;
		if((hits[h1].flag & 0x100) <= 0) continue;
		if((hits[h2].flag & 0x100) <= 0) continue;
		if(primary.find(hits[h1].qname) == primary.end()) continue;
		cnt++;
		redundant[h1] = true;
		redundant[h2] = true;
	}

	printf("filter %d redundant fragments, total %lu frags %lu reads\n", cnt, frgs.size(), hits.size());

	vector<hit> v;
	chain_set s;
	for(int i = 0; i < hits.size(); i++)
	{
		if(redundant[i] == true) continue;
		v.push_back(hits[i]);
		vector<int32_t> chain = hcst.get_chain(i);
		if(chain.size() >= 1) s.add(chain, v.size() - 1);
	}
	hits = v;
	hcst = s;

	return 0;
}

int bundle::build_phase_set(phase_set &ps, splice_graph &gr)
{
	vector<int> fb(hits.size(), -1);
	for(int i = 0; i < frgs.size(); i++)
	{
		if(brdg[i] <= -1) continue;

		int h1 = frgs[i].first;	
		int h2 = frgs[i].second;

		if(brdg[i] == 0)
		{
			fb[h1] = 0;			// paired, to be bridged
			fb[h2] = 0;			// paired, to be bridged
			continue;
		}

		int u1 = gr.locate_vertex(hits[h1].pos);
		int u2 = gr.locate_vertex(hits[h2].rpos - 1);

		if(u1 < 0 || u2 < 0) continue;
		int32_t p1 = gr.get_vertex_info(u1).lpos;
		int32_t p2 = gr.get_vertex_info(u2).rpos;

		vector<int32_t> v1 = hcst.get_chain(h1);
		vector<int32_t> v2 = hcst.get_chain(h2);

		vector<int32_t> xy;

		if(brdg[i] == 1)
		{
			bool b = merge_intron_chains(v1, v2, xy);
			if(b == false) continue;
		}

		if(brdg[i] >= 2)
		{
			vector<int32_t> vv = fcst.get_chain(i);
			xy.insert(xy.end(), v1.begin(), v1.end());
			xy.insert(xy.end(), vv.begin(), vv.end());
			xy.insert(xy.end(), v2.begin(), v2.end());
		}

		xy.insert(xy.begin(), p1);
		xy.insert(xy.end(), p2);

		bool b = check_increasing_sequence(xy);
		if(b == false) continue;

		fb[h1] = 1;			// bridged
		fb[h2] = 1;			// bridged

		ps.add(xy, 1);
	}

	for(int i = 0; i < hits.size(); i++)
	{
		if(fb[i] >= 0) continue;

		int u1 = gr.locate_vertex(hits[i].pos);
		int u2 = gr.locate_vertex(hits[i].rpos - 1);

		if(u1 < 0 || u2 < 0) continue;
		int32_t p1 = gr.get_vertex_info(u1).lpos;
		int32_t p2 = gr.get_vertex_info(u2).rpos;

		vector<int32_t> xy = hcst.get_chain(i);
		xy.insert(xy.begin(), p1);
		xy.insert(xy.end(), p2);

		bool b = check_increasing_sequence(xy);
		if(b == false) continue;

		ps.add(xy, 1);
	}
	return 0;
}

int bundle::update_bridges(const vector<int> &frlist, const vector<int32_t> &chain)
{
	int cnt = 0;
	for(int i = 0; i < frlist.size(); i++)
	{
		assert(chain.size() % 2 == 0);

		int k = frlist[i];

		assert(brdg[k] == 0);
		hit &h1 = hits[frgs[k].first];
		hit &h2 = hits[frgs[k].second];

		vector<int32_t> v1;
		v1.push_back(h1.rpos);
		v1.insert(v1.end(), chain.begin(), chain.end());
		v1.push_back(h2.pos);

		if(h1.rpos < h2.pos && check_increasing_sequence(v1) == false) continue;

		cnt++;

		if(chain.size() <= 0)
		{
			brdg[k] = 1;
		}
		else
		{
			assert(chain.size() >= 2);
			brdg[k] = 2;
			fcst.add(chain, k);
		}

		for(int k = 0; k < v1.size() / 2; k++)
		{
			int32_t p1 = v1[k * 2 + 0];
			int32_t p2 = v1[k * 2 + 1];
			if(p1 >= p2) continue;
			mmap += make_pair(ROI(p1, p2), 1);
		}
	}
	return cnt;
}
