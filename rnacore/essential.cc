#include "util.h"
#include "essential.h"
#include "splice_graph.h"

#include <algorithm>

int transform_vertex_set_map(const set<int> &s, map<int, int> &m)
{
	m.clear();
	if(s.size() <= 0) return 0;

	vector<int> vv(s.begin(), s.end());
	sort(vv.begin(), vv.end());
	for(int i = 0; i < vv.size(); i++)
	{
		m.insert(pair<int, int>(vv[i], i + 1));
	}
	return 0;
}

int build_child_splice_graph(splice_graph &root, splice_graph &gr, map<int, int> &a2b)
{
	gr.clear();
	if(a2b.size() <= 0) return 0;

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
	int32_t x;
	for(int k = 0; k < chain.size() / 2; k++)
	{
		int32_t p = chain[k * 2 + 0];
		int32_t q = chain[k * 2 + 1];
		assert(p < q);
		x += q - p;
	}
	return x;
}

vector<int> project_vector(const vector<int> &v, const map<int, int> &a2b)
{
	vector<int> vv;
	for(int k = 0; k < v.size(); k++)
	{
		map<int, int>::const_iterator it = a2b.find(v[k]);
		if(it == a2b.end()) break;
		vv.push_back(it->second);
	}
	return vv;
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

bool build_path_from_exon_coordinates(splice_graph &gr, const vector<int32_t> &v, vector<int> &vv)
{
	// assume v encodes exon-chain coordinates
	// assume that lindex and rindex are available in gr
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

		// test
		if(a > b)
		{
			printf("k = %d, pp[k] = (%d, %d)\n", k, a, b);
			printf("exon coordinates = ( ");
			printv(v);
			printf(")\n");
			gr.print();
			printf("\n");
		}

		assert(a <= b);
		if(check_continuous_vertices(gr, a, b) == false) return false;
		for(int j = a; j <= b; j++) vv.push_back(j);
	}

	for(int i = 0; i < vv.size() - 1; i++) assert(vv[i] < vv[i + 1]);
	return true;
}

bool build_path_from_intron_coordinates(splice_graph &gr, const vector<int32_t> &v, vector<int> &vv)
{
	// assume v encodes intron chain coordinates
	// assume that lindex and rindex are available in gr
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
	// assume that lindex and rindex are available in gr

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
	// make sure that lindex and rindex are available in gr
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

/*
bool transform_to_paths(splice_graph &gr, PRC &p)
{
	vector<int> v1;
	vector<int> v2;
	bool b1 = build_path_from_exon_coordinates(gr, p.first.vv, v1);
	bool b2 = build_path_from_exon_coordinates(gr, p.second.vv, v2);
	if(b1 == false || b2 == false) return false;
	p.first.vv = v1;
	p.second.vv = v2;
	return true;
}
*/

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
		if(h.isize >= 0) continue;

		// do not use hi; as long as qname, pos and isize are identical
		int k = (h.get_qhash() % max_index + h.pos % max_index + (0 - h.isize) % max_index) % max_index;
		vv[k].push_back(i);
	}

	for(int i = 0; i < hits.size(); i++)
	{
		const hit &h = hits[i];
		if(paired[i] == true) continue;
		if(h.isize <= 0) continue;

		int k = (h.get_qhash() % max_index + h.mpos % max_index + h.isize % max_index) % max_index;
		int x = -1;
		for(int j = 0; j < vv[k].size(); j++)
		{
			int u = vv[k][j];
			const hit &z = hits[u];
			//if(z.hi != h.hi) continue;
			if(paired[u] == true) continue;
			if(z.pos != h.mpos) continue;
			if(z.isize + h.isize != 0) continue;
			//if(z.qhash != h.qhash) continue;
			if(z.qname != h.qname) continue;
			x = vv[k][j];
			break;
		}

		if(x == -1) continue;

		fs.push_back(PI(i, x));

		paired[i] = true;
		paired[x] = true;
	}

	//printf("total hits = %lu, total fragments = %lu\n", hits.size(), fragments.size());
	return 0;
}

int compare_phasing_paths(const vector<int> &ref, const vector<int> &qry)
{
	if(ref.back() < qry.front()) return FALL_RIGHT;
	if(ref.front() > qry.back()) return FALL_LEFT;

	vector<int>::const_iterator r1 = lower_bound(ref.begin(), ref.end(), qry.front());
	vector<int>::const_iterator q1 = lower_bound(qry.begin(), qry.end(), ref.front());
	assert(r1 != ref.end());
	assert(q1 != qry.end());

	int kr1 = r1 - ref.begin();
	int kq1 = q1 - qry.begin();
	assert(kr1 == 0 || kq1 == 0);

	vector<int>::const_iterator q2 = lower_bound(qry.begin(), qry.end(), ref.back());
	vector<int>::const_iterator r2 = lower_bound(ref.begin(), ref.end(), qry.back());
	assert(r2 != ref.end() || q2 != qry.end());

	int kr2 = r2 - ref.begin();
	int kq2 = q2 - qry.begin();

	if(*q1 == ref.front() || *r1 == qry.front())
	{
		if(r2 != ref.end() && q2 != qry.end())
		{
			assert(ref.back() == qry.back());
			assert(kr2 == ref.size() - 1);
			assert(kq2 == qry.size() - 1);
			bool b = identical(ref, kr1, kr2, qry, kq1, kq2);
			if(b == false) return CONFLICTING;
			if(b == true && kr1 == 0 && kq1 == 0) return IDENTICAL;
			if(b == true && kr1 >= 1 && kq1 == 0) return CONTAINED;
			if(b == true && kr1 == 0 && kq1 >= 1) return CONTAINING;
			assert(false);
		}
		else if(r2 != ref.end() && q2 == qry.end())
		{
			bool b = identical(ref, kr1, kr2, qry, kq1, qry.size() - 1);
			if(b == false) return CONFLICTING;
			if(b == true && kq1 == 0) return CONTAINED;
			if(b == true && kq1 >= 1) return EXTEND_LEFT;
			assert(false);
		}
		else if(r2 == ref.end() && q2 != qry.end())
		{
			bool b = identical(ref, kr1, ref.size() - 1, qry, kq1, kq2);
			if(b == false) return CONFLICTING;
			if(b == true && kr1 == 0) return CONTAINING;
			if(b == true && kr1 >= 1) return EXTEND_RIGHT;
		}
	}
	else if(*r1 > qry.front() && r2 == r1 && *r2 > qry.back()) return NESTED;
	else if(*q1 > ref.front() && q2 == q1 && *q2 > ref.back()) return NESTING;
	return CONFLICTING;
}

bool merge_phasing_paths(const vector<int> &ref, const vector<int> &qry, vector<int> &merged)
{
	merged.clear();
	int t = compare_phasing_paths(ref, qry);

	if(t == CONFLICTING) return false;
	if(t == IDENTICAL) merged = ref;
	if(t == FALL_RIGHT) return false;
	if(t == FALL_LEFT) return false;
	if(t == CONTAINED) merged = ref;
	if(t == CONTAINING) merged = qry;
	if(t == NESTED) merged = ref;
	if(t == NESTING) merged = qry;
	if(t == EXTEND_LEFT)
	{
		vector<int>::const_iterator q1 = lower_bound(qry.begin(), qry.end(), ref.front());
		assert(*q1 == ref.front());
		merged.insert(merged.end(), qry.begin(), q1);
		merged.insert(merged.end(), ref.begin(), ref.end());
	}
	if(t == EXTEND_RIGHT)
	{
		vector<int>::const_iterator q2 = lower_bound(qry.begin(), qry.end(), ref.back());
		assert(*q2 == ref.back());
		merged.insert(merged.end(), ref.begin(), ref.end());
		merged.insert(merged.end(), q2 + 1, qry.end());
	}
	return true;
}

bool identical(const vector<int> &x, int x1, int x2, const vector<int> &y, int y1, int y2)
{
	assert(x1 >= 0 && x1 < x.size());
	assert(x2 >= 0 && x2 < x.size());
	assert(y1 >= 0 && y1 < y.size());
	assert(y2 >= 0 && y2 < y.size());

	if(x[x1] != y[y1]) return false;
	if(x[x2] != y[y2]) return false;
	if(x2 - x1 != y2 - y1) return false;

	for(int kx = x1, ky = y1; kx <= x2 && ky <= y2; kx++, ky++)
	{
		if(x[kx] != y[ky]) return false;
	}

	return true;
}
