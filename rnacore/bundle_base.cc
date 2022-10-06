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
#include "essential.h"
#include "constants.h"
#include "util.h"

bundle_base::bundle_base()
{
	tid = -1;
	chrm = "";
	lpos = 1 << 30;
	rpos = 0;
	strand = '.';
	//unbridged = -1;
}

int bundle_base::add_hit_intervals(const hit &ht, bam1_t *b)
{
	add_hit(ht);
	add_intervals(b);
	vector<int32_t> v = ht.extract_splices(b);
	if(v.size() >= 1) hcst.add(v, hits.size() - 1, ht.xs);
	return 0;
}

int bundle_base::add_borrowed_path(const vector<int32_t> &p, double w)
{
	assert(p.size() % 2 == 0);
	for(int k = 0; k < p.size() / 2; k++)
	{
		int p1 = p[k*2+0];
		int p2 = p[k*2+1];
		if(p1 >= 0 && p2 >= 0)
		{
			if(p1 < lpos) lpos = p1;
			if(p2 > rpos) rpos = p2;
			mmap += make_pair(ROI(p1, p2), (int)(w));
		}
		else if(p1 < 0 && p2 < 0)
		{
			vector<int32_t> v;
			v.push_back(0 - p1);
			v.push_back(0 - p2);
			hcst.add(v, -1, strand);
		}
	}
	return 0;
}

int bundle_base::add_hit(const hit &ht)
{
	// store new hit
	hits.push_back(ht);

	// calcuate the boundaries on reference
	if(ht.pos < lpos) lpos = ht.pos;

	int32_t p = ht.rpos;
	// TODO parameter
	if(ht.mpos > ht.rpos && ht.mpos <= ht.rpos + 10000) p = ht.mpos;
	if(p > rpos) rpos = p;

	// set tid
	if(tid == -1) tid = ht.tid;
	assert(tid == ht.tid);

	// set strand
	if(hits.size() <= 1) strand = ht.strand;
	assert(strand == ht.strand);

	return 0;
}

int bundle_base::add_intervals(bam1_t *b)
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
	int n = hits.size();
	int m = frgs.size();
	hits.clear();
	frgs.clear();
	vector<hit>().swap(hits);
	vector<AI3>().swap(frgs);
	hcst.clear();
	fcst.clear();
	mmap.clear();
	imap.clear();
	split_interval_map().swap(mmap);
	split_interval_map().swap(imap);
	return 0;
}

int bundle_base::compute_strand(int libtype)
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

int bundle_base::print(int index)
{
	printf("bundle_base%d: ", index);

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

int bundle_base::build_fragments()
{
	frgs.clear();
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
		if(h.hid < 0) continue;

		// do not use hi; as long as qname, pos and isize are identical
		int k = (h.get_qhash() % max_index + h.pos % max_index + (0 - h.isize) % max_index) % max_index;
		vv[k].push_back(i);
	}

	for(int i = 0; i < hits.size(); i++)
	{
		const hit &h = hits[i];
		if(h.hid < 0) continue;
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
		frgs.push_back(AI3({i, x, 0}));
		paired[i] = true;
		paired[x] = true;
	}

	//printf("total hits = %lu, total fragments = %lu\n", hits.size(), frgs.size());
	return 0;
}

int bundle_base::count_unbridged()
{
	int unbridged = 0;
	for(int i = 0; i < frgs.size(); i++)
	{
		// only group unbridged fragments
		if(frgs[i][2] >= 1) continue;
		if(frgs[i][2] <= -1) continue;
		unbridged++;
	}
	return unbridged;
}

int bundle_base::build_phase_set(phase_set &ps, splice_graph &gr)
{
	vector<int> fb(hits.size(), -1);
	for(int i = 0; i < frgs.size(); i++)
	{
		if(frgs[i][2] <= -1) continue;

		int h1 = frgs[i][0];	
		int h2 = frgs[i][1];

		assert(hits[h1].hid >= 0);
		assert(hits[h2].hid >= 0);

		if(frgs[i][2] == 0)
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

		vector<int32_t> v1 = hcst.get(h1).first;
		vector<int32_t> v2 = hcst.get(h2).first;

		vector<int32_t> xy;

		if(frgs[i][2] == 1)
		{
			bool b = merge_intron_chains(v1, v2, xy);
			if(b == false) continue;
		}

		if(frgs[i][2] >= 2)
		{
			vector<int32_t> vv = fcst.get(i).first;
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
		if(hits[i].hid < 0) continue;

		int u1 = gr.locate_vertex(hits[i].pos);
		int u2 = gr.locate_vertex(hits[i].rpos - 1);

		if(u1 < 0 || u2 < 0) continue;
		int32_t p1 = gr.get_vertex_info(u1).lpos;
		int32_t p2 = gr.get_vertex_info(u2).rpos;

		vector<int32_t> xy = hcst.get(i).first;
		xy.insert(xy.begin(), p1);
		xy.insert(xy.end(), p2);

		bool b = check_increasing_sequence(xy);
		if(b == false) continue;

		ps.add(xy, 1);
	}
	return 0;
}

int bundle_base::update_bridges(const vector<int> &frlist, const vector<int32_t> &chain)
{
	int cnt = 0;
	for(int i = 0; i < frlist.size(); i++)
	{
		assert(chain.size() % 2 == 0);

		int k = frlist[i];
		assert(frgs[k][2] == 0);

		hit &h1 = hits[frgs[k][0]];
		hit &h2 = hits[frgs[k][1]];
		
		assert(h1.hid >= 0);
		assert(h2.hid >= 0);

		// possibly revise pos/rpos
		/*
		if(chain.size() >= 1 && h1.rpos > chain.front())
		{
			printf("retract small region 1: %d-%d\n", chain.front(), h1.rpos);
			mmap -= make_pair(ROI(chain.front(), h1.rpos), 1);
			h1.rpos = chain.front();
		}
		if(chain.size() >= 1 && h2.pos < chain.back())
		{
			printf("retract small region 2: %d-%d\n", h2.pos, chain.back());
			mmap -= make_pair(ROI(h2.pos, chain.back()), 1);
			h2.pos = chain.back();
		}
		*/

		vector<int32_t> v1;
		v1.push_back(h1.rpos);

		v1.insert(v1.end(), chain.begin(), chain.end());
		v1.push_back(h2.pos);

		if(h1.rpos < h2.pos && check_increasing_sequence(v1) == false) continue;

		cnt++;

		if(chain.size() <= 0)
		{
			frgs[k][2] = 1;
		}
		else
		{
			assert(chain.size() >= 2);
			frgs[k][2] = 2;
			if(h1.xs == h2.xs) fcst.add(chain, k, h1.xs);
			else fcst.add(chain, k, '.');
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

int bundle_base::eliminate_bridge(int k)
{
	assert(k >= 0 && k < frgs.size());
	assert(frgs[k][2] >= 1);

	hit &h1 = hits[frgs[k][0]];
	hit &h2 = hits[frgs[k][1]];
	assert(h1.hid >= 0);
	assert(h2.hid >= 0);

	vector<int32_t> chain = fcst.get(k).first;

	vector<int32_t> v1;
	v1.push_back(h1.rpos);
	v1.insert(v1.end(), chain.begin(), chain.end());
	v1.push_back(h2.pos);

	for(int i = 0; i < v1.size() / 2; i++)
	{
		int32_t p1 = v1[i * 2 + 0];
		int32_t p2 = v1[i * 2 + 1];
		if(p1 >= p2) continue;
		mmap += make_pair(ROI(p1, p2), -1);
	}

	frgs[k][2] = -1;
	fcst.remove(k);

	return 0;
}

int bundle_base::eliminate_hit(int k)
{
	assert(k >= 0 && k < hits.size());

	hit &h1 = hits[k];
	assert(h1.hid >= 0);

	vector<int32_t> chain = hcst.get(k).first;

	vector<int32_t> v1;
	v1.push_back(h1.pos);
	v1.insert(v1.end(), chain.begin(), chain.end());
	v1.push_back(h1.rpos);

	for(int i = 0; i < v1.size() / 2; i++)
	{
		int32_t p1 = v1[i * 2 + 0];
		int32_t p2 = v1[i * 2 + 1];
		if(p1 >= p2) continue;
		mmap += make_pair(ROI(p1, p2), -1);
	}

	h1.hid = -1;
	hcst.remove(k);

	return 0;
}

int bundle_base::filter_secondary_hits()
{
	set<string> primary;
	for(int i = 0; i < frgs.size(); i++)
	{
		int h1 = frgs[i][0];
		int h2 = frgs[i][1];
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
		int h1 = frgs[i][0];
		int h2 = frgs[i][1];
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
		vector<int32_t> chain = hcst.get(i).first;
		if(chain.size() >= 1) s.add(chain, v.size() - 1, hits[i].xs);
	}
	hits = v;
	hcst = s;

	return 0;
}

int bundle_base::filter_multialigned_hits()
{
	set<string> bridged;
	set<string> primary;
	for(int i = 0; i < frgs.size(); i++)
	{
		if(frgs[i][2] <= 0) continue;
		int h1 = frgs[i][0];
		int h2 = frgs[i][1];
		assert(hits[h1].qname == hits[h2].qname);
		bridged.insert(hits[h1].qname);
		if((hits[h1].flag & 0x100) <= 0 && (hits[h2].flag & 0x100) <= 0) primary.insert(hits[h1].qname);
	}

	int cnt1 = 0;

	// remove unbridged pairs
	for(int i = 0; i < frgs.size(); i++)
	{
		int h1 = frgs[i][0];
		int h2 = frgs[i][1];
		if(frgs[i][2] >= 1) continue;
		if(primary.find(hits[h1].qname) == primary.end()) continue;
		eliminate_hit(h1);
		eliminate_hit(h2);
		frgs[i][2] = -1;
		cnt1++;
	}

	// remove bridged but secondary pairs
	for(int i = 0; i < frgs.size(); i++)
	{
		int h1 = frgs[i][0];
		int h2 = frgs[i][1];
		if(frgs[i][2] <= 0) continue;
		if((hits[h1].flag & 0x100) <= 0) continue;
		if((hits[h2].flag & 0x100) <= 0) continue;
		if(primary.find(hits[h1].qname) == primary.end()) continue;
		eliminate_bridge(i);
		eliminate_hit(h1);
		eliminate_hit(h2);
		cnt1++;
	}

	// mark unpaired hits
	vector<bool> paired(hits.size(), false);
	for(int i = 0; i < frgs.size(); i++)
	{
		int h1 = frgs[i][0];
		int h2 = frgs[i][1];
		paired[h1] = true;
		paired[h2] = true;
	}

	// remove unpaired hits
	int cnt2 = 0;
	for(int i = 0; i < hits.size(); i++)
	{
		if(paired[i] == true) continue;
		if(bridged.find(hits[i].qname) == bridged.end()) continue;
		eliminate_hit(i);
		cnt2++;
	}

	//printf("filter %d / %d multialigned fragments / hits, total %lu frags %lu reads\n", cnt1, cnt2, frgs.size(), hits.size());
	return 0;
}
