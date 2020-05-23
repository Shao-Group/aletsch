/*
Part of meta-scallop
(c) 2020 by Mingfu Shao, The Pennsylvania State University
See LICENSE for licensing.
*/

#include "util.h"
#include "essential.h"
#include "splice_graph.h"

#include <algorithm>
#include <cstdio>
#include <cstring>

int build_child_splice_graph(splice_graph &root, splice_graph &gr, map<int, int> &a2b)
{
	gr.clear();
	if(a2b.size() <= 0) return 0;

	/*
	printf("constructing child graph: ");
	for(auto &z : a2b) printf(" %d -> %d, ", z.first, z.second);
	printf("\n");
	*/

	vector<int> vv = get_keys(a2b);
	sort(vv.begin(), vv.end());

	gr.chrm = root.chrm;
	gr.strand = root.strand;

	int32_t lpos = root.get_vertex_info(vv.front()).lpos;
	int32_t rpos = root.get_vertex_info(vv.back()).rpos;

	// vertices
	gr.add_vertex();
	vertex_info vi0;
	vi0.lpos = lpos;
	vi0.rpos = lpos;
	gr.set_vertex_weight(0, 0);
	gr.set_vertex_info(0, vi0);

	for(int i = 0; i < vv.size(); i++)
	{
		int k = vv[i];
		gr.add_vertex();
		gr.set_vertex_weight(i + 1, root.get_vertex_weight(k));
		gr.set_vertex_info(i + 1, root.get_vertex_info(k));
	}

	gr.add_vertex();
	vertex_info vin;
	vin.lpos = rpos;
	vin.rpos = rpos;
	gr.set_vertex_weight(vv.size() + 1, 0);
	gr.set_vertex_info(vv.size() + 1, vin);

	// edges
	edge_iterator it1, it2;
	PEEI pei;
	for(pei = root.out_edges(0), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		int t = (*it1)->target();

		//printf("there exists edge from 0 to %d in root\n", t);
		if(a2b.find(t) == a2b.end()) continue;
		int y = a2b[t];
		//printf("y = a2b[%d] = %d\n", t, y);

		edge_descriptor e = gr.add_edge(0, y);
		gr.set_edge_weight(e, root.get_edge_weight(*it1));
		gr.set_edge_info(e, root.get_edge_info(*it1));
	}

	int n = root.num_vertices() - 1;
	for(int i = 0; i < vv.size(); i++)
	{
		int s = vv[i];
		assert(s != 0 && s != n);
		assert(a2b.find(s) != a2b.end());
		int x = a2b[s];

		pei = root.out_edges(s); 
		for(it1 = pei.first; it1 != pei.second; it1++)
		{
			int t = (*it1)->target();
			assert(t == n || a2b.find(t) != a2b.end());
			int y = ((t == n) ? gr.num_vertices() - 1 : a2b[t]);

			//printf("there exists edge from %d to %d in root => (%d, %d)\n", s, t, x, y);

			edge_descriptor e = gr.add_edge(x, y);
			gr.set_edge_weight(e, root.get_edge_weight(*it1));
			gr.set_edge_info(e, root.get_edge_info(*it1));
		}
	}
	return 0;
}

int32_t get_total_length_of_introns(const vector<int32_t> &chain)
{
	assert(chain.size() % 2 == 0);
	int32_t x = 0;
	for(int k = 0; k < chain.size() / 2; k++)
	{
		int32_t p = chain[k * 2 + 0];
		int32_t q = chain[k * 2 + 1];
		assert(p < q);
		x += q - p;
	}
	return x;
}

int build_exon_coordinates_from_path(splice_graph &gr, const vector<int> &v, vector<int32_t> &vv)
{
	vv.clear();
	if(v.size() <= 0) return 0;

	int n = gr.num_vertices() - 1;
	int32_t pre = -99999;

	if(v.front() == 0) vv.push_back(-1);
	if(v.front() == 0) vv.push_back(-1);

	for(int i = 0; i < v.size(); i++)
	{
		int p = v[i];

		if(p == 0) continue;
		if(p == n) continue;

		int32_t pp = gr.get_vertex_info(p).lpos;
		int32_t qq = gr.get_vertex_info(p).rpos;

		if(pp == pre) 
		{
			pre = qq;
			continue;
		}

		if(pre >= 0) vv.push_back(pre);
		vv.push_back(pp);

		pre = qq;
	}

	if(pre >= 0) vv.push_back(pre);
	if(v.back() == n) vv.push_back(-2);
	if(v.back() == n) vv.push_back(-2);

	return 0;
}

int build_intron_coordinates_from_path(splice_graph &gr, const vector<int> &v, vector<int32_t> &vv)
{
	vv.clear();
	for(int i = 0; i < v.size() - 1; i++)
	{
		int32_t pp = gr.get_vertex_info(v[i + 0]).rpos;
		int32_t qq = gr.get_vertex_info(v[i + 1]).lpos;

		assert(pp <= qq);
		if(pp == qq) continue;
		vv.push_back(pp);
		vv.push_back(qq);
	}
	return 0;
}

int check_strand_from_intron_coordinates(splice_graph &gr, const vector<int32_t> &v)
{
	// assume v encodes intron chain coordinates
	assert(v.size() % 2 == 0);
	if(v.size() <= 0) return 0;

	bool b1 = false;
	bool b2 = false;
	int n = v.size() / 2;
	for(int k = 0; k < n; k++)
	{
		int32_t p = v[2 * k + 0];
		int32_t q = v[2 * k + 1];
		assert(p >= 0 && q >= 0);
		if(p >= q) return -1;
		if(gr.rindex.find(p) == gr.rindex.end()) return -1;
		if(gr.lindex.find(q) == gr.lindex.end()) return -1;
		int kp = gr.rindex[p];
		int kq = gr.lindex[q];

		PEB pe = gr.edge(kp, kq);
		if(pe.second == false) return -1;

		int strand = gr.get_edge_info(pe.first).strand;
		if(strand == 1) b1 = true;
		if(strand == 2) b2 = true;
	}

	if(b1 == true && b2 == true) return -1;
	if(b1 == true) return 1;
	if(b2 == true) return 2;
	return 0;
}

bool build_path_from_exon_coordinates(splice_graph &gr, const vector<int32_t> &v, vector<int> &vv)
{
	// assume v encodes exon-chain coordinates
	vv.clear();
	assert(v.size() % 2 == 0);
	if(v.size() <= 0) return true;

	int n = v.size() / 2;
	vector<PI> pp(n);
	for(int k = 0; k < n; k++)
	{
		int32_t p = v[2 * k + 0];
		int32_t q = v[2 * k + 1];
		assert(p >= 0 && q >= 0);
		assert(p <= q);

		if(gr.lindex.find(p) == gr.lindex.end()) return false;
		if(gr.rindex.find(q) == gr.rindex.end()) return false;
		int kp = gr.lindex[p];
		int kq = gr.rindex[q];
		pp[k].first = kp;
		pp[k].second = kq;
	}

	for(int k = 0; k < n; k++)
	{
		int a = pp[k].first;
		int b = pp[k].second;

		if(a > b)
		{
			printf("TEST: a = %d, b = %d\n", a, b);
			gr.print();
			printf("v = ");
			printv(v);
			printf("\n list = ");
			for(int i = 0; i < n; i++) printf("(%d, %d), ", pp[i].first, pp[i].second);
			printf("\n");
			fflush(stdout);
		}

		if(a > b) return false;
		//assert(a <= b);

		if(check_continuous_vertices(gr, a, b) == false) return false;
		for(int j = a; j <= b; j++) vv.push_back(j);
	}

	for(int i = 0; i < vv.size() - 1; i++) assert(vv[i] < vv[i + 1]);
	return true;
}

bool build_path_from_intron_coordinates(splice_graph &gr, const vector<int32_t> &v, vector<int> &vv)
{
	// assume v encodes intron chain coordinates
	vv.clear();
	assert(v.size() % 2 == 0);
	if(v.size() <= 0) return true;

	int n = v.size() / 2;
	vector<PI> pp(n);
	for(int k = 0; k < n; k++)
	{
		int32_t p = v[2 * k + 0];
		int32_t q = v[2 * k + 1];
		assert(p >= 0 && q >= 0);
		assert(p <= q);

		if(gr.rindex.find(p) == gr.rindex.end()) return false;
		if(gr.lindex.find(q) == gr.lindex.end()) return false;
		int kp = gr.rindex[p];
		int kq = gr.lindex[q];
		pp[k].first = kp;
		pp[k].second = kq;
	}

	vv.push_back(pp.front().first);
	for(int k = 0; k < n - 1; k++)
	{
		int a = pp[k + 0].second;
		int b = pp[k + 1].first;
		assert(a <= b);
		if(check_continuous_vertices(gr, a, b) == false) return false;
		for(int j = a; j <= b; j++) vv.push_back(j);
	}
	vv.push_back(pp.back().second);

	return true;
}

bool build_path_from_mixed_coordinates(splice_graph &gr, const vector<int32_t> &v, vector<int> &vv)
{
	// assume v[1..n-1] encodes intron-chain coordinates
	vv.clear();
	assert(v.size() % 2 == 0);
	if(v.size() <= 0) return true;

	int u1 = gr.locate_vertex(v.front());
	int u2 = gr.locate_vertex(v.back() - 1);

	if(u1 < 0 || u2 < 0) return false;

	if(v.size() == 2)
	{
		for(int k = u1; k <= u2; k++) vv.push_back(k);
		return true;
	}

	vector<int> uu;
	vector<int32_t> u(v.begin() + 1, v.end() - 1);
	bool b = build_path_from_intron_coordinates(gr, u, uu);
	
	if(b == false) return false;

	for(int i = u1; i < uu.front(); i++) vv.push_back(i);
	vv.insert(vv.end(), uu.begin(), uu.end());
	for(int i = uu.back() + 1; i <= u2; i++) vv.push_back(i);

	return true;
}

bool check_continuous_vertices(splice_graph &gr, int x, int y)
{
	if(x >= y) return true;
	for(int i = x; i < y; i++)
	{
		PEB p = gr.edge(i, i + 1);
		if(p.second == false) return false;
		if(gr.get_vertex_info(i).rpos != gr.get_vertex_info(i + 1).lpos) return false;
	}
	return true;
}

bool check_valid_path(splice_graph &gr, const vector<int> &vv)
{
	int n = gr.num_vertices() - 1;
	for(int k = 0; k < vv.size() - 1; k++)
	{
		if(vv[k + 0] < 0 || vv[k + 0] > n) return false;
		if(vv[k + 1] < 0 || vv[k + 1] > n) return false;
		PEB p = gr.edge(vv[k], vv[k + 1]);
		if(p.second == false) return false;
	}
	return true;
}

bool align_hit_to_splice_graph(const hit &h, splice_graph &gr, vector<int> &vv)
{
	vv.clear();
	vector<int64_t> v;
	h.get_aligned_intervals(v);

	if(v.size() == 0) return false;

	vector<int32_t> u;
	for(int i = 0; i < v.size(); i++)
	{
		u.push_back(high32(v[i]));
		u.push_back(low32(v[i]));
	}

	bool b = build_path_from_mixed_coordinates(gr, u, vv);
	return b;
}

int build_paired_reads(const vector<hit> &hits, vector<PI> &fs)
{
	vector<bool> paired(hits.size(), false);

	fs.clear();
	if(hits.size() == 0) return 0;

	int max_index = hits.size() + 1;
	if(max_index > 1000000) max_index = 1000000;

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
		fs.push_back(PI(i, x));

		paired[i] = true;
		paired[x] = true;
	}

	//printf("total hits = %lu, total fragments = %lu\n", hits.size(), fs.size());
	return 0;
}

bool merge_intron_chains(const vector<int32_t> &x, const vector<int32_t> &y, vector<int32_t> &xy)
{
	xy.clear();
	if(x.size() >= 1 && y.size() >= 1 && x.front() > y.front()) return false;
	bool b = merge_two_sorted_sequences<int32_t>(x, y, xy);
	if(b == false) return false;
	int d = x.size() + y.size() - xy.size();
	if(d % 2 != 0) return false;
	return true;
}

bool consistent_intron_chains(const vector<int32_t> &x, const vector<int32_t> &y)
{
	vector<int32_t> v;
	return merge_intron_chains(x, y, v);
}

int add_cigar_match(bam1_t &b1t, int32_t p1, int32_t p2)
{
	assert(p1 < p2);
	uint32_t c1 = p2 - p1;
	c1 = c1 << 4;
	c1 += BAM_CMATCH;
	memcpy(b1t.data + b1t.l_data, &c1, 4);
	b1t.l_data += 4;
	b1t.core.n_cigar++;
	return 0;
}

int add_cigar_skip(bam1_t &b1t, int32_t p1, int32_t p2)
{
	assert(p1 < p2);
	uint32_t c1 = p2 - p1;
	c1 = c1 << 4;
	c1 += BAM_CREF_SKIP;
	memcpy(b1t.data + b1t.l_data, &c1, 4);
	b1t.l_data += 4;
	b1t.core.n_cigar++;
	return 0;
}

bool build_bam1_t(bam1_t &b1t, const hit_core &h, const vector<int32_t> &chain)
{
	b1t.core = h;
	b1t.core.pos = h.pos;
	b1t.core.bin = 0;
	b1t.core.n_cigar = 0;
	b1t.core.l_qseq = 0;
	b1t.core.isize = h.isize;
	b1t.core.flag = h.flag;

	b1t.m_data = b1t.core.l_qname + 4 * (chain.size() + 1) + 7 * 3;
	b1t.data = new uint8_t[b1t.m_data];

	// copy qname
	b1t.l_data = 0;
	assert(h.qname.size() + b1t.core.l_extranul + 1 == b1t.core.l_qname);
	memcpy(b1t.data, h.qname.c_str(), h.qname.size());
	b1t.l_data += h.qname.length();
	for(int i = 0; i <= b1t.core.l_extranul; i++)
	{
		b1t.data[b1t.l_data] = 0;
		b1t.l_data++;
	}
	assert(b1t.l_data == b1t.core.l_qname);

	vector<int32_t> z;
	z.push_back(h.pos);
	z.insert(z.end(), chain.begin(), chain.end());
	z.push_back(h.rpos);

	// CIGAR
	for(int i = 0; i < z.size() - 1; i++)
	{
		int32_t x1 = z[i + 0];
		int32_t x2 = z[i + 1];

		if(x1 >= x2)
		{
			printf("fail to build bam1 for %s\n", h.qname.c_str());
			return false;
		}

		if(i % 2 == 0) add_cigar_match(b1t, x1, x2);
		else add_cigar_skip(b1t, x1, x2);
	}

	if(h.xs !='.') bam_aux_append(&(b1t), "XS", 'A', 1, (uint8_t*)(&h.xs));
	if(h.hi != -1) bam_aux_append(&(b1t), "HI", 'i', 4, (uint8_t*)(&h.hi));
	if(h.nh != -1) bam_aux_append(&(b1t), "NH", 'i', 4, (uint8_t*)(&h.nh));

	return true;
}

bool build_bam1_t(bam1_t &b1t, const hit_core &h1, const hit_core &h2, const vector<int32_t> &chain)
{
	b1t.core = h1;
	b1t.core.pos = h1.pos;
	b1t.core.bin = 0;
	b1t.core.n_cigar = 0;
	b1t.core.l_qseq = 0;
	b1t.core.isize = h2.rpos - h1.pos;
	b1t.core.flag = h1.flag;
	//b1t.core.flag -= b1t.core.flag & (0x1);

	b1t.m_data = b1t.core.l_qname + 4 * (chain.size() + 1) + 7 * 3;
	b1t.data = new uint8_t[b1t.m_data];

	// copy qname
	b1t.l_data = 0;
	assert(h1.qname == h2.qname);
	assert(h1.qname.size() + b1t.core.l_extranul + 1 == b1t.core.l_qname);
	memcpy(b1t.data, h1.qname.c_str(), h1.qname.size());
	b1t.l_data += h1.qname.length();
	for(int i = 0; i <= b1t.core.l_extranul; i++)
	{
		b1t.data[b1t.l_data] = 0;
		b1t.l_data++;
	}
	assert(b1t.l_data == b1t.core.l_qname);

	vector<int32_t> z;
	z.push_back(h1.pos);
	z.insert(z.end(), chain.begin(), chain.end());
	z.push_back(h2.rpos);

	// CIGAR
	for(int i = 0; i < z.size() - 1; i++)
	{
		int32_t x1 = z[i + 0];
		int32_t x2 = z[i + 1];

		if(x1 >= x2)
		{
			printf("fail to build bam1 for %s\n", h1.qname.c_str());
			return false;
		}

		//assert(x1 <= x2);
		if(i % 2 == 0) add_cigar_match(b1t, x1, x2);
		else add_cigar_skip(b1t, x1, x2);
	}

	char c = h1.xs;
	if(c == '.' && h2.xs != '.') c = h2.xs;

	if(c != '.') bam_aux_append(&(b1t), "XS", 'A', 1, (uint8_t*)(&c));
	if(h1.hi != -1 && h1.hi == h2.hi) bam_aux_append(&(b1t), "HI", 'i', 4, (uint8_t*)(&h1.hi));
	if(h1.nh != -1 && h1.nh == h2.nh) bam_aux_append(&(b1t), "NH", 'i', 4, (uint8_t*)(&h1.nh));

	return true;
}

int write_bridged_pereads_cluster(BGZF *fout, const pereads_cluster &pc, const vector<int32_t> &whole)
{
	assert(pc.hits1.size() == pc.hits2.size());
	if(pc.hits1.size() == 0) return 0;

	for(int i = 0; i < pc.hits1.size(); i++)
	{
		const hit &h1 = pc.hits1[i];
		const hit &h2 = pc.hits2[i];

		bam1_t b1t;

		/*
		printf("write bam1 (%d, %d) -- (%d, %d), whole = ( ", h1.pos, h1.rpos, h2.pos, h2.rpos);
		printv(whole);
		printf("\n");
		*/

		bool b = build_bam1_t(b1t, h1, h2, whole);
		if(b == true) bam_write1(fout, &(b1t));
		// TODO, handle false
		assert(b1t.data != NULL);
		delete b1t.data;
	}
	return 0;
}

int write_unbridged_pereads_cluster(BGZF *fout, const pereads_cluster &pc)
{
	assert(pc.hits1.size() == pc.hits2.size());
	if(pc.hits1.size() == 0) return 0;

	for(int i = 0; i < pc.hits1.size(); i++)
	{
		const hit &h1 = pc.hits1[i];

		bam1_t b1t;
		bool b = build_bam1_t(b1t, h1, pc.chain1);
		if(b == true) bam_write1(fout, &(b1t));
		assert(b1t.data != NULL);
		delete b1t.data;
	}

	for(int i = 0; i < pc.hits2.size(); i++)
	{
		const hit &h2 = pc.hits2[i];

		bam1_t b1t;
		bool b = build_bam1_t(b1t, h2, pc.chain2);
		if(b == true) bam_write1(fout, &(b1t));
		assert(b1t.data != NULL);
		delete b1t.data;
	}
	return 0;
}
