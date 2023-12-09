/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#include "splice_graph.h"
#include "util.h"
#include "config.h"
#include "interval_map.h"
#include <sstream>
#include <fstream>
#include <cfloat>
#include <cmath>
#include <algorithm>
#include <cstdlib>

using namespace std;

splice_graph::splice_graph()
{}

splice_graph::splice_graph(const splice_graph &gr)
{
	chrm = gr.chrm;
	gid = gr.gid;
	strand = gr.strand;
    reads = gr.reads;
    subgraph = gr.subgraph;
	lindex = gr.lindex;
	rindex = gr.rindex;

	MEE x2y;
	MEE y2x;
	copy(gr, x2y, y2x);
}

int splice_graph::copy(const splice_graph &gr, MEE &x2y, MEE &y2x)
{
	clear();
	for(int i = 0; i < gr.num_vertices(); i++)
	{
		add_vertex();
		set_vertex_weight(i, gr.get_vertex_weight(i));
		set_vertex_info(i, gr.get_vertex_info(i));
	}

	PEEI p = gr.edges();
	for(edge_iterator it = p.first; it != p.second; it++)
	{
		edge_descriptor e = add_edge((*it)->source(), (*it)->target());
		set_edge_weight(e, gr.get_edge_weight(*it));
		set_edge_info(e, gr.get_edge_info(*it));

		assert(e != NULL);
		assert(ewrt.find(e) != ewrt.end());
		assert(einf.find(e) != einf.end());
		assert(x2y.find(*it) == x2y.end());
		assert(y2x.find(e) == y2x.end());

		x2y.insert(PEE(*it, e));
		y2x.insert(PEE(e, *it));
	}

	return 0;
}

int splice_graph::clear()
{
	directed_graph::clear();
	vwrt.clear();
	vinf.clear();
	ewrt.clear();
	einf.clear();
	lindex.clear();
	rindex.clear();
	return 0;
}

splice_graph::~splice_graph()
{}

double splice_graph::get_vertex_weight(int v) const
{
	assert(v >= 0 && v < vwrt.size());
	return vwrt[v];
}

vertex_info & splice_graph::get_editable_vertex_info(int v)
{
	assert(v >= 0 && v < vinf.size());
	return vinf[v];
}

const vertex_info & splice_graph::get_vertex_info(int v) const
{
	assert(v >= 0 && v < vinf.size());
	return vinf[v];
}

double splice_graph::get_edge_weight(edge_base *e) const
{
	MED::const_iterator it = ewrt.find(e);
	assert(it != ewrt.end());
	return it->second;
}

edge_info & splice_graph::get_editable_edge_info(edge_base *e)
{
	MEIF::iterator it = einf.find(e);
	assert(it != einf.end());
	return it->second;
}

const edge_info & splice_graph::get_edge_info(edge_base *e) const
{
	MEIF::const_iterator it = einf.find(e);
	assert(it != einf.end());
	return it->second;
}

int splice_graph::set_vertex_weight(int v, double w) 
{
	assert(v >= 0 && v < vv.size());
	if(vwrt.size() != vv.size()) vwrt.resize(vv.size());
	vwrt[v] = w;
	return 0;
}

int splice_graph::set_vertex_info(int v, const vertex_info &vi) 
{
	assert(v >= 0 && v < vv.size());
	if(vinf.size() != vv.size()) vinf.resize(vv.size());
	vinf[v] = vi;
	return 0;
}

int splice_graph::set_edge_weight(edge_base* e, double w) 
{
	if(ewrt.find(e) != ewrt.end()) ewrt[e] = w;
	else ewrt.insert(PED(e, w));
	return 0;
}

int splice_graph::set_edge_info(edge_base* e, const edge_info &ei) 
{
	if(einf.find(e) != einf.end()) einf[e] = ei;
	else einf.insert(PEIF(e, ei));
	return 0;
}

MED splice_graph::get_edge_weights() const
{
	return ewrt;
}

vector<double> splice_graph::get_vertex_weights() const
{
	return vwrt;
}

int splice_graph::set_edge_weights(const MED &med)
{
	ewrt = med;
	return 0;
}

int splice_graph::set_vertex_weights(const vector<double> &v)
{
	vwrt = v;
	return 0;
}

double splice_graph::get_out_weights(int v)
{
	PEEI pei;
	edge_iterator it1, it2;
	double ww = 0;
	for(pei = out_edges(v), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		double w = get_edge_weight(*it1);
		ww += w;
	}
	return ww;
}

double splice_graph::get_in_weights(int v)
{
	PEEI pei;
	edge_iterator it1, it2;
	double ww = 0;
	for(pei = in_edges(v), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		double w = get_edge_weight(*it1);
		ww += w;
	}
	return ww;
}

double splice_graph::get_max_out_weight(int v)
{
	edge_iterator it1, it2;
	PEEI pei;
	double ww = 0;
	for(pei = out_edges(v), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		double w = get_edge_weight(*it1);
		if(w < ww) continue;
		ww = w;
	}
	return ww;
}

double splice_graph::get_max_in_weight(int v)
{
	edge_iterator it1, it2;
	PEEI pei;
	double ww = 0;
	for(pei = in_edges(v), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		double w = get_edge_weight(*it1);
		if(w < ww) continue;
		ww = w;
	}
	return ww;
}

edge_descriptor splice_graph::max_out_edge(int v)
{
	edge_iterator it1, it2;
	edge_descriptor ee = null_edge;
	PEEI pei;
	double ww = 0;
	for(pei = out_edges(v), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		double w = get_edge_weight(*it1);
		if(w < ww) continue;
		ee = (*it1);
		ww = w;
	}
	return ee;
}

edge_descriptor splice_graph::max_in_edge(int v)
{
	PEEI pei;
	edge_iterator it1, it2;
	edge_descriptor ee = null_edge;
	double ww = 0;
	for(pei = in_edges(v), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		double w = get_edge_weight(*it1);
		if(w < ww) continue;
		ee = (*it1);
		ww = w;
	}
	return ee;
}

edge_descriptor splice_graph::min_out_edge(int v)
{
	edge_iterator it1, it2;
	edge_descriptor ee = null_edge;
	PEEI pei;
	double ww = DBL_MAX;
	for(pei = out_edges(v), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		double w = get_edge_weight(*it1);
		if(w > ww) continue;
		ee = (*it1);
		ww = w;
	}
	return ee;
}

edge_descriptor splice_graph::min_in_edge(int v)
{
	PEEI pei;
	edge_iterator it1, it2;
	edge_descriptor ee = null_edge;
	double ww = DBL_MAX;
	for(pei = in_edges(v), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		double w = get_edge_weight(*it1);
		if(w > ww) continue;
		ee = (*it1);
		ww = w;
	}
	return ee;
}

int splice_graph::count_junctions()
{
	int cnt = 0;
	edge_iterator it1, it2;
	PEEI pei;
	for(pei = edges(), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		int s = (*it1)->source();
		int t = (*it1)->target();
		if(s == 0) continue;
		if(t == num_vertices() - 1) continue;
		int32_t p1 = get_vertex_info(s).rpos;
		int32_t p2 = get_vertex_info(t).lpos;
		if(p1 >= p2) continue;
		double w = get_edge_weight(*it1);
		cnt += (int)(w);
	}
	return cnt;
}

edge_descriptor splice_graph::compute_maximum_edge_w()
{
	edge_iterator it1, it2;
	double max_weight = -1;
	edge_descriptor max_edge = null_edge;
	PEEI pei;

	for(pei = edges(), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		double w = get_edge_weight(*it1);
		if(w < max_weight) continue;
		max_weight = w;
		max_edge = *it1;
	}
	return max_edge;
}

int splice_graph::build(const string &file)
{
	ifstream fin(file.c_str());
	if(fin.fail()) 
	{
		printf("open file %s error\n", file.c_str());
		return 0;
	}

	char line[10240];
	// get the number of vertices
	fin.getline(line, 10240, '\n');	
	int n = atoi(line);

	for(int i = 0; i < n; i++)
	{
		char name[10240];
		double weight;
		vertex_info vi;
		fin.getline(line, 10240, '\n');	
		stringstream sstr(line);
		sstr>>name>>weight>>vi.length;

		add_vertex();
		set_vertex_weight(i, weight);
		set_vertex_info(i, vi);
	}

	while(fin.getline(line, 10240, '\n'))
	{
		int x, y;
		double weight;
		edge_info ei;
		stringstream sstr(line);
		sstr>>x>>y>>weight>>ei.length;

		assert(x != y);
		assert(x >= 0 && x < num_vertices());
		assert(y >= 0 && y < num_vertices());

		edge_descriptor p = add_edge(x, y);
		set_edge_weight(p, weight);
		set_edge_info(p, ei);
	}

	fin.close();
	return 0;
}

int splice_graph::write(ostream &os, const vector<int> &v, double w, const string &name)
{
	if(v.size() <= 0) return 0;
	int n = num_vertices() - 1;
	vector<int32_t> vv;
	int32_t pre = -9999;

	if(v.front() == 0) vv.push_back(-1);
	if(v.front() == 0) vv.push_back(-1);

	for(int i = 0; i < v.size(); i++)
	{
		int p = v[i];

		if(p == 0) continue;
		if(p == n) continue;

		int32_t pp = get_vertex_info(p).lpos;
		int32_t qq = get_vertex_info(p).rpos;

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

	if(vv.size() <= 1) return 0;

	os << name.c_str() << " " << vv.size();
	for(int k = 0; k < vv.size(); k++) os << " " << vv[k];
	os << " " << w << " 1" << endl;
	return 0;
}

int splice_graph::write(ostream &os)
{
	os<<fixed;
	os.precision(2);
	int n = num_vertices() - 1;
	
	os << "# " << gid.c_str() << " " << chrm.c_str() << " " << strand << endl;
	for(int i = 1; i < n; i++)
	{
		double weight = get_vertex_weight(i);
		vertex_info vi = get_vertex_info(i);
		if(vi.lpos >= vi.rpos) continue;
		os<<"region " << vi.lpos << " " << vi.rpos << " " << weight << endl;
	}

	PEEI pei = out_edges(0);
	for(edge_iterator it = pei.first; it != pei.second; it++)
	{
		int s = (*it)->source(); 
		int t = (*it)->target();
		assert(s == 0 && t > s);
		if(t == n) continue;
		double w = get_edge_weight(*it);
		int32_t p = get_vertex_info(t).lpos;
		os << "sbound " << p << " " << w << " 1" << endl;
	}

	pei = in_edges(n);
	for(edge_iterator it = pei.first; it != pei.second; it++)
	{
		int s = (*it)->source(); 
		int t = (*it)->target();
		assert(t == n);
		assert(s < t);
		if(s == 0) continue;
		double w = get_edge_weight(*it);
		int32_t p = get_vertex_info(s).rpos;
		os << "tbound " << p << " " << w << " 1" << endl;
	}

	pei = edges();
	for(edge_iterator it = pei.first; it != pei.second; it++)
	{
		int s = (*it)->source(); 
		int t = (*it)->target();
		assert(s < t);
		if(s == 0) continue;
		if(t == n) continue;
		double w = get_edge_weight(*it);
		int32_t p1 = get_vertex_info(s).rpos;
		int32_t p2 = get_vertex_info(t).lpos;
		if(p1 >= p2) continue;
		os << "junction " << p1 << " " << p2 << " " << w << " 1" << endl;
	}
	return 0;
}

int splice_graph::simulate(int nv, int ne, int mw)
{
	clear();
	for(int i = 0; i < nv; i++)
	{
		add_vertex();
		vinf.push_back(vertex_info());
		vwrt.push_back(0);
	}

	while(true)
	{
		int s = rand() % nv;
		if(s == nv - 1) continue;
		int t = rand() % ((nv - s - 1) / 2 + 1) + s + 1;
		assert(s >= 0 && s < nv);
		assert(t > s  && t < nv);
		if(s == 0 && t == nv - 1) continue;
		int f = rand() % (mw - 10) + 10;
		PEB p = edge(s, t);
		if(p.second == true) continue;

		edge_descriptor e = add_edge(s, t);
		ewrt.insert(PED(e, f));
		einf.insert(PEIF(e, edge_info()));
		if(ewrt.size() >= ne) break;
	}

	assert(in_degree(0) == 0);
	assert(out_degree(num_vertices() - 1) == 0);

	MED med;
	while(true)
	{
		VE v;
		int w = (int)(compute_maximum_path_w(v));
		if(w <= 0) break;
		for(int i = 0; i < v.size(); i++)
		{
			ewrt[v[i]] -= w;
			if(med.find(v[i]) == med.end()) med.insert(PED(v[i], w));
			else med[v[i]] += w;
		}
	}

	for(MED::iterator it = ewrt.begin(); it != ewrt.end(); it++)
	{
		if(med.find(it->first) == med.end()) remove_edge(it->first);
	}

	ewrt = med;
	einf.clear();
	for(MED::iterator it = ewrt.begin(); it != ewrt.end(); it++)
	{
		einf.insert(PEIF(it->first, edge_info()));
	}

	edge_iterator it1, it2;
	PEEI pei;
	for(int i = 0; i < num_vertices(); i++)
	{
		int wx = 0;
		for(pei = in_edges(i), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
		{
			wx += (int)(ewrt[*it1]);
		}
		int wy = 0;
		for(pei = out_edges(i), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
		{
			wy += (int)(ewrt[*it1]);
		}

		if(i == 0) assert(wx == 0);
		else if(i == num_vertices() - 1) wy = 0;
		else assert(wx == wy);

		if(i == 0) vwrt[i] = wy;
		else vwrt[i] = wx;
	}

	assert(vwrt[0] == vwrt[num_vertices() - 1]);

	int vv = 0;
	for(int i = 0; i < num_vertices(); i++)
	{
		if(degree(i) >= 1) vv++;
	}
	int delta = num_edges() - vv + 2;
	printf("simulate %d vertices, %lu edges, %.0lf max-flow, delta = %d\n", vv, num_edges(), vwrt[0], delta);

	return 0;
}

int splice_graph::compute_decomp_paths()
{
	int n = 0;
	for(int i = 0; i < num_vertices(); i++)
	{
		if(degree(i) == 0) continue;
		n++;
	}
	if(n == 0) return 0;
	int m = num_edges();
	return (m - n + 2);
}

long splice_graph::compute_num_paths()
{
	long max = 9999999999;
	vector<long> table;
	int n = num_vertices();
	table.resize(n, 0);
	table[0] = 1;
	for(int i = 1; i < n; i++)
	{
		edge_iterator it1, it2;
		PEEI pei;
		for(pei = in_edges(i), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
		{
			int s = (*it1)->source();
			int t = (*it1)->target();
			assert(t == i);
			//assert(s < i);
			table[t] += table[s];
			if(table[t] >= max) return max;
		}
	}
	
	return table[n - 1];
}

long splice_graph::compute_num_paths(int a, int b)
{
	long max = 9999999999;
	vector<long> table;
	int n = (b - a + 1);
	table.resize(n, 0);
	table[0] = 1;
	for(int i = a + 1; i <= b; i++)
	{
		edge_iterator it1, it2;
		PEEI pei;
		for(pei = in_edges(i), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
		{
			int s = (*it1)->source();
			int t = (*it1)->target();
			assert(t == i);
			if(s < a) continue;
			table[t] += table[s - a];
			if(table[t] >= max) return max;
		}
	}
	return table[n - 1];
}

bool splice_graph::check_fully_connected()
{
	assert(num_vertices() >= 2);
	if(num_vertices() <= 2) return true;

	vector<int> s;
	vector<int> t;
	bfs(0, s);
	bfs_reverse(num_vertices() - 1, t);
	
	if(s.size() != num_vertices()) return false;
	if(t.size() != num_vertices()) return false;
	return true;
}

int splice_graph::compute_independent_subgraphs()
{
	int num = 0;
	edge_iterator it1, it2;
	PEEI pei;
	set<int> ss;
	ss.insert(num_vertices() - 1);
	for(pei = in_edges(num_vertices() - 1), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		int x = (*it1)->source();
		ss.insert(x);
	}
	
	for(pei = out_edges(0), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		int x = (*it1)->target();
		while(true)
		{
			if(ss.find(x) != ss.end()) num++;
			if(ss.find(x) != ss.end()) break;
			x = compute_out_equivalent_vertex(x);
			if(x == -1) break;
		}
	}
	return num;
}

int splice_graph::bfs_w(int s, double w, vector<int> &v, VE &b)
{
	v.assign(num_vertices(), -1);
	b.assign(num_vertices(), null_edge);
	vector<bool> closed;
	closed.resize(num_vertices(), false);
	vector<int> open;
	open.push_back(s);
	v[s] = 0;
	int p = 0;

	while(p < open.size())
	{
		int x = open[p];
		assert(v[x] >= 0);

		p++;
		if(closed[x] == true) continue;
		closed[x] = true;

		edge_iterator it1, it2;
		PEEI pei;
		for(pei = out_edges(x), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
		{
			int y = (*it1)->target();
			double ww = get_edge_weight(*it1);
			//if(ww < w - SMIN) continue;
			if(ww < w) continue;
			if(v[y] == -1) 
			{
				v[y] = 1 + v[x];
				b[y] = (*it1);
			}
			assert(v[y] <= 1 + v[x]);
			if(closed[y] == true) continue;
			open.push_back(y);
		}
	}
	return 0;
}

int splice_graph::compute_closest_path(int s, vector<double> &d)
{
	vector<int> b;
	return compute_closest_path(s, d, b);
}

int splice_graph::compute_closest_path_reverse(int t, vector<double> &d)
{
	vector<int> b;
	return compute_closest_path_reverse(t, d, b);
}

int splice_graph::compute_closest_path(int s, vector<double> &d, vector<int> &b)
{
	int n = num_vertices() - s;
	d.assign(n, -1);
	b.assign(n, -1);
	d[0] = 0;
	for(int i = s + 1; i < num_vertices(); i++)
	{
		double w0 = -1;
		int b0 = -1;
		edge_iterator it1, it2;
		PEEI pei;
		for(pei = in_edges(i), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
		{
			int ss = (*it1)->source();
			if(ss < s) continue;
			int k = ss - s;
			double ww = get_edge_weight(*it1);
			if(b0 == -1 || d[k] + ww < w0)
			{
				w0 = d[k] + ww;
				b0 = k;
			}
		}
		d[i - s] = w0;
		b[i - s] = b0;
	}
	return 0;
}

int splice_graph::compute_closest_path_reverse(int t, vector<double> &d, vector<int> &b)
{
	d.assign(t + 1, -1);
	b.assign(t + 1, -1);
	d[t] = 0;
	for(int i = t - 1; i >= 0; i--)
	{
		double w0 = -1;
		int b0 = -1;
		edge_iterator it1, it2;
		PEEI pei;
		for(pei = out_edges(i), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
		{
			int tt = (*it1)->target();
			if(tt > t) continue;
			double ww = get_edge_weight(*it1);
			if(b0 == -1 || d[tt] + ww < w0)
			{
				w0 = d[tt] + ww;
				b0 = tt;
			}
		}
		d[i] = w0;
		b[i] = b0;
	}
	return 0;
}

int splice_graph::compute_shortest_path_w(int s, int t, double w)
{
	vector<int> v;
	VE b;
	bfs_w(s, w, v, b);
	return v[t];
}

int splice_graph::compute_shortest_path_w(int s, int t, double w, VE &p)
{
	vector<int> v;
	VE b;
	bfs_w(s, w, v, b);
	if(v[t] == -1) return -1;

	p.clear();
	int x = t;
	while(x != s)
	{
		assert(b[x] != null_edge);
		p.push_back(b[x]);
		x = b[x]->source();
	}
	reverse(p.begin(), p.end());
	return v[t];
}

double splice_graph::compute_maximum_path_w(VE &p)
{
	return compute_maximum_st_path_w(p, 0, num_vertices() - 1);
}

double splice_graph::compute_maximum_st_path_w(VE &p, int ss, int tt)
{
	p.clear();
	vector<double> table;		// dynamic programming table
	VE back;					// backtrace edge pointers
	table.resize(num_vertices(), -1);
	back.resize(num_vertices(), null_edge);

	vector<int> tp = topological_sort();
	int n = num_vertices();
	assert(tp.size() == n);
	//assert(tp[0] == 0);
	//assert(tp[n - 1] == n - 1);

	int ssi = -1;
	int tti = -1;
	for(int i = 0; i < tp.size(); i++)
	{
		if(tp[i] == ss) ssi = i;
		if(tp[i] == tt) tti = i;
	}
	assert(ssi != -1);
	assert(tti != -1);

	table[ss] = DBL_MAX;
	for(int ii = ssi + 1; ii <= tti; ii++)
	{
		int i = tp[ii];
		if(degree(i) == 0) continue;

		double max_abd = 0;
		edge_descriptor max_edge = null_edge;
		edge_iterator it1, it2;
		PEEI pei;
		for(pei = in_edges(i), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
		{
			int s = (*it1)->source();
			int t = (*it1)->target();
			assert(t == i);
			if(table[s] <= -1) continue;
			double xw = get_edge_weight(*it1);
			double ww = xw < table[s] ? xw : table[s];
			if(ww >= max_abd)
			{
				max_abd = ww;
				max_edge = *it1;
			}
		}

		if(max_edge == null_edge) continue;

		back[i] = max_edge;
		table[i] = max_abd;
	}

	int x = tt;
	while(true)
	{
		edge_descriptor e = back[x]; 
		if(e == null_edge) break;
		p.push_back(e);
		x = e->source();
	}
	reverse(p.begin(), p.end());

	return table[tt];
}

bool splice_graph::compute_optimal_path(VE &p)
{
	p.clear();
	vector<double> table;		// dynamic programming table
	VE back;					// backtrace edge pointers
	table.resize(num_vertices(), 0);
	back.resize(num_vertices(), null_edge);
	table[0] = 0;

	vector<int> tp = topological_sort();
	int n = num_vertices();
	assert(tp.size() == n);
	assert(tp[0] == 0);
	assert(tp[n - 1] == n - 1);

	edge_iterator it1, it2;
	PEEI pei;
	for(int ii = 1; ii < n; ii++)
	{
		int i = tp[ii];
		if(degree(i) == 0) continue;

		double sum = 0;
		for(pei = in_edges(i), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
		{
			sum += get_edge_weight(*it1);
		}

		double req = DBL_MAX;
		edge_descriptor ee = null_edge;
		for(pei = in_edges(i), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
		{
			int s = (*it1)->source();
			int t = (*it1)->target();
			double xw = get_edge_weight(*it1);
			double ww = sum - xw + table[s];
			if(ww < req)
			{
				req = ww;
				ee = *it1;
			}
		}
		assert(ee != null_edge);

		back[i] = ee;
		table[i] = req;
	}

	edge_descriptor opt = null_edge;
	for(pei = in_edges(n - 1), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		int s = (*it1)->source();
		double w = get_edge_weight(*it1);
		if(w < table[s] + 1) continue;
		opt = *it1;
		break;
	}

	if(opt == null_edge) return false;

	edge_descriptor e = opt;
	while(e != null_edge)
	{
		p.push_back(e);
		int s = e->source();
		e = back[s];
	}
	reverse(p.begin(), p.end());

	return true;
}

double splice_graph::compute_minimum_weight(const VE &p)
{
	double min = DBL_MAX;
	for(int i = 0; i < p.size(); i++)
	{
		double w = get_edge_weight(p[i]);
		if(w < min) min = w;
	}
	return min;
}

int splice_graph::locate(int v)
{
	if(v == 0) return 1;
	if(v == num_vertices() - 1) return 2;
	edge_iterator it1, it2;
	PEEI pei;
	bool b1 = false, b2 = false;
	for(pei = in_edges(v), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		if((*it1)->source() == 0) b1 = true;
	}
	for(pei = out_edges(v), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		if((*it1)->target() == num_vertices() - 1) b2 = true;
	}

	if(b1 == true && b2 == true) return 3;
	if(b1 == true) return 4;
	if(b2 == true) return 5;
	return 0;
}

int splice_graph::round_weights()
{
	MED m = ewrt;
	for(MED::iterator it = m.begin(); it != m.end(); it++)
	{
		it->second = 0.0;
	}

	while(true)
	{
		edge_descriptor e = compute_maximum_edge_w();
		if(e == null_edge) break;

		double w0 = get_edge_weight(e);
		if(w0 <= 0) break;

		VE v1, v2;
		double w1 = w0;
		double w2 = w0;

		if(e->source() != 0) w1 = compute_maximum_st_path_w(v1, 0, e->source());
		if(e->target() != num_vertices() - 1) w2 = compute_maximum_st_path_w(v2, e->target(), num_vertices() - 1);

		assert(w1 <= w0);
		assert(w2 <= w0);

		double w = (w1 < w2) ? w1 : w2;
		double ww = ceil(w);
		if(ww <= 0) ww = 1;

		VE v = v1;
		v.push_back(e);
		v.insert(v.end(), v2.begin(), v2.end());
		
		for(int i = 0; i < v.size(); i++)
		{
			m[v[i]] += ww;
			ewrt[v[i]] -= ww;
			if(ewrt[v[i]] <= 0) ewrt[v[i]] = 0;
		}
	}

	ewrt = m;
	vwrt.assign(num_vertices(), 0);
	edge_iterator it1, it2;
	PEEI pei;
	for(pei = out_edges(0), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		double w = ewrt[*it1];
		vwrt[0] += w;
	}

	for(int i = 1; i < num_vertices(); i++)
	{
		for(pei = in_edges(i), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
		{
			double w = ewrt[*it1];
			vwrt[i] += w;
		}
	}

	return 0;
}

double splice_graph::compute_average_edge_weight()
{
	edge_iterator it1, it2;
	PEEI pei;
	int cnt = 0;
	double sum = 0;
	for(pei = edges(), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		int s = (*it1)->source();
		int t = (*it1)->target();
		if(s == 0) continue;
		if(t == num_vertices() - 1) continue;
		cnt++;
		sum += get_edge_weight(*it1);
	}
	if(cnt >= 1) sum = sum / cnt;
	return sum;
}

double splice_graph::compute_average_vertex_weight()
{
	double sum = 0;
	int cnt = num_vertices() - 2;
	for(int i = 1; i < num_vertices() - 1; i++)
	{
		sum += get_vertex_weight(i);
	}
	if(cnt >= 1) sum = sum / cnt;
	return sum;
}

int splice_graph::build_vertex_index()
{
	lindex.clear();
	rindex.clear();
	int n = num_vertices() - 1;
	for(int i = 0; i <= n; i++)
	{
		const vertex_info &v = get_vertex_info(i);
		if(i != 0) lindex.insert(pair<int32_t, int>(v.lpos, i));
		if(i != n) rindex.insert(pair<int32_t, int>(v.rpos, i));
	}
	return 0;
}

int splice_graph::determine_position_right_type(int32_t p)
{
	if(rindex.find(p) == rindex.end()) return -1;
	int x = rindex[p];
	int n = num_vertices() - 1;

	if(edge(x, n).second == true) return END_BOUNDARY;

	bool splice = false;
	bool adjacent = false;
	PEEI pei = out_edges(x);
	vertex_info vx = get_vertex_info(x);
	for(edge_iterator it = pei.first; it != pei.second; it++)
	{
		edge_descriptor e = *it;
		vertex_info vi = get_vertex_info(e->target());
		if(vi.lpos == vx.rpos) adjacent = true;
		if(vi.lpos >  vx.rpos) splice = true;
	}

	if(splice == true && adjacent == true) return RIGHT_MIXED;
	if(splice == true) return RIGHT_SPLICE;
	if(adjacent == true && determine_position_left_type(p) == LEFT_MIXED) return LEFT_MIXED;

	return MIDDLE_CUT;
}

int splice_graph::determine_position_left_type(int32_t p)
{
	if(lindex.find(p) == lindex.end()) return -1;
	int x = lindex[p];

	if(edge(0, x).second == true) return START_BOUNDARY;

	bool splice = false;
	bool adjacent = false;
	PEEI pei = in_edges(x);
	vertex_info vx = get_vertex_info(x);
	for(edge_iterator it = pei.first; it != pei.second; it++)
	{
		edge_descriptor e = *it;
		vertex_info vi = get_vertex_info(e->source());
		if(vi.rpos == vx.lpos) adjacent = true;
		if(vi.rpos <  vx.lpos) splice = true;
	}

	if(splice == true && adjacent == true) return LEFT_MIXED;
	if(splice == true) return LEFT_SPLICE;
	if(adjacent == true && determine_position_right_type(p) == RIGHT_MIXED) return RIGHT_MIXED;

	return MIDDLE_CUT;
}

int splice_graph::locate_lbound(int32_t p)
{
	if(lindex.find(p) == lindex.end()) return -1;
	else return lindex[p];
}

int splice_graph::locate_rbound(int32_t p)
{
	if(rindex.find(p) == rindex.end()) return -1;
	else return rindex[p];
}

int splice_graph::locate_vertex(int32_t p)
{
	int m = locate_vertex(p, 1, num_vertices() - 1);
	assert(m >= 1 && m <= num_vertices() - 1);
	if(m == num_vertices()) return -1;
	const vertex_info &v = get_vertex_info(m);
	if(p >= v.lpos && p < v.rpos) return m;
	else return -1;
}

int splice_graph::locate_vertex0(int32_t p, int a, int b)
{
	if(a >= b) return -1;
	int m = (a + b) / 2;
	assert(m >= 0 && m < num_vertices());
	const vertex_info &v = get_vertex_info(m);
	if(p >= v.lpos && p < v.rpos) return m;
	if(p < v.lpos) return locate_vertex(p, a, m);
	return locate_vertex(p, m + 1, b);
}

int splice_graph::locate_vertex(int32_t p, int a, int b)
{
	// return the leftmost vertex v such that p < v.rpos
	// when this function is called, vertex b (exclusive)
	// won't be checked
	if(a >= b) return b;
	int m = (a + b) / 2;
	assert(m >= a && m < b);
	assert(m >= 0 && m < num_vertices());
	const vertex_info &v = get_vertex_info(m);
	if(p >= v.lpos && p < v.rpos) return m;
	if(p < v.lpos) return locate_vertex(p, a, m);
	return locate_vertex(p, m + 1, b);
}


int splice_graph::draw(const string &file, const MIS &mis, const MES &mes, double len, const vector<int> &tp, bool footer)
{
	return directed_graph::draw(file, mis, mes, len, tp, footer);
}

int splice_graph::draw(const string &file, const MIS &mis, const MES &mes, double len, bool footer)
{
	return directed_graph::draw(file, mis, mes, len, footer);
}

int splice_graph::draw(const string &file, bool footer)
{
	MIS mis;
	char buf[10240];

	for(int i = 0; i < num_vertices(); i++)
	{
		double w = get_vertex_weight(i);
		vertex_info vi = get_vertex_info(i);
		int ll = vi.lpos % 100000;
		int rr = vi.rpos % 100000;
		sprintf(buf, "1:%.1lf:%d:%d:%d", w, rr - ll, ll, rr);
		mis.insert(PIS(i, buf));
	}

	MES mes;
	edge_iterator it1, it2;
	PEEI pei;
	for(pei = edges(), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		double w = get_edge_weight(*it1);
		sprintf(buf, "%.1lf:%d", w, get_edge_info(*it1).strand);
		mes.insert(PES(*it1, buf));
	}
	draw(file, mis, mes, 1.5, footer);
	return 0;
}

int splice_graph::print_nontrivial_vertices()
{
	int k = 0;
	for(int i = 1; i < num_vertices() - 1; i++)
	{
		if(in_degree(i) <= 1) continue;
		if(out_degree(i) <= 1) continue;
		printf("nontrivial vertex %d length = %d, weight = %.2lf\n", k++, get_vertex_info(i).length, get_vertex_weight(i));
	}
	return 0;
}

int splice_graph::print()
{
	for(int i = 0; i < num_vertices(); i++)
	{
		vertex_info vi = get_vertex_info(i);
		printf("vertex %d, range = [%d, %d), length = %d, weight = %.2lf\n", i, vi.lpos, vi.rpos, vi.rpos - vi.lpos, get_vertex_weight(i));
	}
	print_weights();
	return 0;
}

int splice_graph::print_weights()
{
	edge_iterator it1, it2;
	PEEI pei;
	for(pei = edges(), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		edge_descriptor e = (*it1);
		int s = e->source();
		int t = e->target();
		int32_t p1 = get_vertex_info(s).rpos;
		int32_t p2 = get_vertex_info(t).lpos;
		double w1 = get_edge_weight(e);
		double w2 = get_edge_info(e).weight;
		int strand = get_edge_info(e).strand;
		printf("edge (%d, %d) pos = %d-%d length = %d weight = (%.2lf, %.2lf) strand = %d\n", s, t, p1, p2, p2 - p1 + 1, w1, w2, strand);
	}
	return 0;
}

int splice_graph::print_junction_supports()
{
    
	edge_iterator it1, it2;
	PEEI pei;
	for(pei = edges(), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
        edge_descriptor e = (*it1);
        int s = e->source();
        int t = e->target();
        int32_t p1 = get_vertex_info(s).rpos;
        int32_t p2 = get_vertex_info(t).lpos;
        double w1 = get_edge_weight(e);

        edge_info ei = get_edge_info(e);
        double w2 = ei.weight;
        if(ei.count >= 0)
        {
            printf("edge (%d, %d) pos = %d-%d length = %d weight = (%.2lf, %.2lf) strand = %d count = %d\n", s, t, p1, p2, p2 - p1 + 1, w1, w2, ei.strand, ei.count);
            printf("support set: ");
            for (auto sp : ei.samples)
                printf("%d(%.2lf), ", sp, ei.spAbd[sp]);
            printf("\n");
        }
    }
    return 0;
}

int32_t splice_graph::get_total_length_of_vertices(const vector<int>& v) const
{
	if(v.size() == 0) return 0;
	int32_t flen = 0;
	for(int i = 0; i < v.size(); i++)
	{
		assert(v[i] >= 0 && v[i] < num_vertices());
		const vertex_info &x = get_vertex_info(v[i]);
		flen += x.rpos - x.lpos;
	}
	return flen;
}

int splice_graph::extend_strands()
{
	PEEI pei = edges();
	for(edge_iterator it = pei.first; it != pei.second; it++)
	{
		edge_descriptor e = *it;
		int sd = get_edge_info(e).strand;
		int s = e->source();
		int t = e->target();
		int32_t p1 = get_vertex_info(s).rpos;
		int32_t p2 = get_vertex_info(t).lpos;
		if(p1 >= p2) continue;
		if(s + 2 != t) continue;

		double we = get_edge_weight(e);
		double wv = get_vertex_weight(s + 1);
		if(we <= wv) continue;
		if(get_vertex_info(s + 1).lpos != p1) continue;
		if(get_vertex_info(s + 1).rpos != p2) continue;
		PEB e1 = edge(s, s + 1);
		PEB e2 = edge(s + 1, t);
		if(e1.second == true) 
		{
			edge_info ei = get_edge_info(e1.first);
			if(ei.strand == 0) ei.strand = sd;
			set_edge_info(e1.first, ei);
		}
		if(e2.second == true) 
		{
			edge_info ei = get_edge_info(e2.first);
			if(ei.strand == 0) ei.strand = sd;
			set_edge_info(e2.first, ei);
		}
	}
	return 0;
}

bool splice_graph::mixed_strand_vertex(int i)
{
	vector<int> v = get_strand_degree(i);
	int p = v[1] + v[4];
	int q = v[2] + v[5];
	if(p >= 1 && q >= 1) return true;
	else return false;
}

vector<int> splice_graph::get_strand_degree(int i)
{
	vector<int> vs(6, 0);
	PEEI pei = in_edges(i);
	for(edge_iterator it = pei.first; it != pei.second; it++)
	{
		edge_descriptor e = *it;
		int s = get_edge_info(e).strand;
		assert(s >= 0 && s <= 2);
		vs[s]++;
	}

	pei = out_edges(i);
	for(edge_iterator it = pei.first; it != pei.second; it++)
	{
		edge_descriptor e = *it;
		int s = get_edge_info(e).strand;
		assert(s >= 0 && s <= 2);
		vs[s + 3]++;
	}

	return vs;
}

bool splice_graph::mixed_strand_graph()
{
	bool p = false;
	bool q = false;
	PEEI pei = edges();
	for(edge_iterator it = pei.first; it != pei.second; it++)
	{
		edge_descriptor e = *it;
		int s = get_edge_info(e).strand;
		if(s == 1) p = true;
		if(s == 2) q = true;
		if(p && q) return true;
	}
	return false;
}

int splice_graph::stat_strandness()
{
	for(int i = 1; i < num_vertices() - 1; i++)
	{
		int32_t p1 = get_vertex_info(i).lpos;
		int32_t p2 = get_vertex_info(i).rpos;
		PEEI pi = in_edges(i);
		PEEI po = out_edges(i);
		int xi[3] = {0, 0, 0};
		int xo[3] = {0, 0, 0};
		double wi[3] = {0, 0, 0};
		double wo[3] = {0, 0, 0};
		for(edge_iterator it = pi.first; it != pi.second; it++)
		{
			edge_descriptor e = (*it);
			edge_info ei = get_edge_info(e);
			double w = get_edge_weight(e);
			xi[ei.strand]++;
			wi[ei.strand] += w;
		}
		for(edge_iterator it = po.first; it != po.second; it++)
		{
			edge_descriptor e = (*it);
			edge_info ei = get_edge_info(e);
			double w = get_edge_weight(e);
			xo[ei.strand]++;
			wo[ei.strand] += w;
		}
		
		if( (xi[1] <= 0 || xo[2] <= 0) && (xi[2] <= 0 || xo[1] <= 0) ) continue;

		printf("vertex %d: %s:%d-%d, in-degree = %d / %d / %d, out-degree = %d / %d / %d, in-weight = %.0lf / %.0lf / %.0lf, out-weight = %.0lf / %.0lf / %.0lf\n",
				i, chrm.c_str(), p1, p2, xi[0], xi[1], xi[2], xo[0], xo[1], xo[2], wi[0], wi[1], wi[2], wo[0], wo[1], wo[2]);
	}
	return 0;
}

int splice_graph::print_vertex(int i)
{
	int32_t p1 = get_vertex_info(i).lpos;
	int32_t p2 = get_vertex_info(i).rpos;
	PEEI pi = in_edges(i);
	PEEI po = out_edges(i);
	int xi[3] = {0, 0, 0};
	int xo[3] = {0, 0, 0};
	double wi[3] = {0, 0, 0};
	double wo[3] = {0, 0, 0};
	for(edge_iterator it = pi.first; it != pi.second; it++)
	{
		edge_descriptor e = (*it);
		edge_info ei = get_edge_info(e);
		double w = get_edge_weight(e);
		xi[ei.strand]++;
		wi[ei.strand] += w;
	}
	for(edge_iterator it = po.first; it != po.second; it++)
	{
		edge_descriptor e = (*it);
		edge_info ei = get_edge_info(e);
		double w = get_edge_weight(e);
		xo[ei.strand]++;
		wo[ei.strand] += w;
	}

	printf("vertex %d: %s:%d-%d, in-degree = %d / %d / %d, out-degree = %d / %d / %d, in-weight = %.0lf / %.0lf / %.0lf, out-weight = %.0lf / %.0lf / %.0lf\n",
			i, chrm.c_str(), p1, p2, xi[0], xi[1], xi[2], xo[0], xo[1], xo[2], wi[0], wi[1], wi[2], wo[0], wo[1], wo[2]);
	return 0;
}

