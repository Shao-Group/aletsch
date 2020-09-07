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
}

int bundle_base::add_hit_intervals(const hit &ht, bam1_t *b)
{
	add_hit(ht);
	add_intervals(b);
	vector<int32_t> v = ht.extract_splices(b);
	if(v.size() >= 1) hcst.add(v, hits.size() - 1, ht.xs);
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
	hits.clear();
	frgs.clear();
	grps.clear();
	hcst.clear();
	fcst.clear();
	pcst.clear();
	qcst.clear();
	gcst.clear();
	mmap.clear();
	imap.clear();
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

	printf("tid = %d, range = %s:%d-%d, orient = %c, #hits = %lu, +/-/. = %d / %d / %d, frgs = %lu, grps = %lu, mmap = %lu, imap = %lu.",
			tid, chrm.c_str(), lpos, rpos, strand, hits.size(), np, nq, n0, frgs.size(), grps.size(), mmap.size(), imap.size());
	hcst.print_size();
	fcst.print_size();
	pcst.print_size();
	qcst.print_size();
	gcst.print_size();
	printf("\n");

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

int bundle_base::group_fragments(splice_graph &gr, int max_gap)
{
	vector<vector<int>> groups;
	group_fragments(gr, groups);
	for(int k = 0; k < groups.size(); k++)
	{
		group_fragments(groups[k], max_gap);
	}
	return 0;
}

int bundle_base::group_fragments(splice_graph &gr, vector<vector<int>> &groups)
{
	typedef pair< vector<int>, vector<int> > PVV;
	map<PVV, int> findex;
	groups.clear();

	for(int i = 0; i < frgs.size(); i++)
	{
		frgs[i][2] = -1;

		int h1 = frgs[i][0];
		int h2 = frgs[i][1];

		assert(hits[h1].hid >= 0);
		assert(hits[h2].hid >= 0);

		if(hits[h1].pos > hits[h2].pos) continue;
		if(hits[h1].rpos > hits[h2].rpos) continue;

		vector<int> v1;
		vector<int> v2;
		const vector<int32_t> &chain1 = hcst.get_chain(h1);
		const vector<int32_t> &chain2 = hcst.get_chain(h2);

		// TODO, don't group with a graph
		bool b1 = align_hit_to_splice_graph(hits[h1], chain1, gr, v1);
		bool b2 = align_hit_to_splice_graph(hits[h2], chain2, gr, v2);

		if(b1 == false || b2 == false)  continue;
		if(v1.size() == 0 || v2.size() == 0) continue;

		PVV pvv(v1, v2);
		if(findex.find(pvv) == findex.end())
		{
			vector<int> v;
			v.push_back(i);
			findex.insert(pair<PVV, int>(pvv, groups.size()));
			groups.push_back(v);
		}
		else
		{
			int k = findex[pvv];
			groups[k].push_back(i);
		}
	}

	//for(int k = 0; k < groups.size(); k++) printf("group %d contains %lu frags\n", k, groups[k].size());

	return 0;
}

int bundle_base::group_fragments(const vector<int> &fs, int max_gap)
{
	vector<vector<int32_t>> vv;
	for(int i = 0; i < fs.size(); i++)
	{
		int h1 = frgs[fs[i]][0];
		int h2 = frgs[fs[i]][1];

		vector<int32_t> v;
		v.push_back(hits[h1].pos);
		v.push_back(hits[h1].rpos);
		v.push_back(hits[h2].pos);
		v.push_back(hits[h2].rpos);
		v.push_back(i);
		vv.push_back(v);
	}

	vector<vector<int>> zz = partition(vv, 0, max_gap);

	for(int i = 0; i < zz.size(); i++)
	{
		if(zz[i].size() <= 0) continue;

		int g = grps.size();

		int h1 = frgs[fs[zz[i][0]]][0];
		int h2 = frgs[fs[zz[i][0]]][1];
		assert(hits[h1].rpos <= hits[h2].rpos);
		assert(hits[h1].pos <= hits[h2].pos);

		AI6 pc;
		const vector<int32_t> &c1 = hcst.get_chain(h1);
		const vector<int32_t> &c2 = hcst.get_chain(h2);
		int s1 = hcst.get_strand(h1);
		int s2 = hcst.get_strand(h2);

		if(c1.size() >= 1) pcst.add(c1, g, s1);
		if(c2.size() >= 1) qcst.add(c2, g, s2);

		vector<int32_t> bounds(4, 0);
		bounds[0] = hits[h1].pos;
		bounds[1] = hits[h1].rpos;
		bounds[2] = hits[h2].pos;
		bounds[3] = hits[h2].rpos;

		for(int k = 0; k < zz[i].size(); k++)
		{
			int f = fs[zz[i][k]];
			h1 = frgs[f][0];
			h2 = frgs[f][1];

			pc[0] += hits[h1].pos  - bounds[0];
			pc[1] += hits[h1].rpos - bounds[1];
			pc[2] += hits[h2].pos  - bounds[2];
			pc[3] += hits[h2].rpos - bounds[3];

			frgs[f][2] = g;
		}

		pc[0] = pc[0] / pc[4] + bounds[0];  
		pc[1] = pc[1] / pc[4] + bounds[1];
		pc[2] = pc[2] / pc[4] + bounds[2]; 
		pc[3] = pc[3] / pc[4] + bounds[3];
		pc[4] = zz[i].size();
		pc[5] = 0;

		grps.push_back(std::move(pc));
	}
	return 0;
}

vector<vector<int>> bundle_base::partition(vector<vector<int32_t>> &fs, int r, int max_partition_gap)
{
	vector<vector<int>> vv;
	if(fs.size() == 0) return vv;

	if(r >= 4)
	{
		vector<int> v;
		for(int k = 0; k < fs.size(); k++) v.push_back(fs[k][4]);
		vv.push_back(v);
		return vv;
	}

	if(r == 0) sort(fs.begin(), fs.end(), compare_rank0);	
	if(r == 1) sort(fs.begin(), fs.end(), compare_rank1);	
	if(r == 2) sort(fs.begin(), fs.end(), compare_rank2);	
	if(r == 3) sort(fs.begin(), fs.end(), compare_rank3);	

	int pre = 0;
	for(int k = 1; k <= fs.size(); k++)
	{
		if(k < fs.size()) assert(fs[k][r] >= fs[k - 1][r]);
		if(k < fs.size() && fs[k][r] - fs[k - 1][r] <= max_partition_gap) continue;

		vector< vector<int32_t> > fs1;
		for(int i = pre; i < k; i++) fs1.push_back(fs[i]);

		vector< vector<int> > vv1 = partition(fs1, r + 1, max_partition_gap);
		vv.insert(vv.end(), vv1.begin(), vv1.end());

		pre = k;
	}
	return vv;
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

bool compare_rank0(const vector<int32_t> &x, const vector<int32_t> &y) { return x[0] < y[0]; }
bool compare_rank1(const vector<int32_t> &x, const vector<int32_t> &y) { return x[1] < y[1]; }
bool compare_rank2(const vector<int32_t> &x, const vector<int32_t> &y) { return x[2] < y[2]; }
bool compare_rank3(const vector<int32_t> &x, const vector<int32_t> &y) { return x[3] < y[3]; }
