#include "predecompose.h"
#include "config.h"

#include <cstdio>
#include <iostream>
#include <climits>
#include <cfloat>
#include <algorithm>

predecompose::predecompose(const splice_graph &g)
	: gr(g)
{
	gr.get_edge_indices(i2e, e2i);
	init_super_edges();
	init_vertex_map();
	init_inner_weights();
	init_nonzeroset();
}

int predecompose::assemble()
{
	while(true)
	{	
		int c = resolve_trivial_vertex();
		if(c <= 0) break;
		else continue;
	}
	collect_paths();
	return 0;
}

int predecompose::resolve_trivial_vertex()
{
	int cnt = 0;
	vector<int> vv(nonzeroset.begin(), nonzeroset.end());
	for(int k = 0; k < vv.size(); k++)
	{
		int i = vv[k];
		assert(gr.in_degree(i) >= 1);
		assert(gr.out_degree(i) >= 1);
		if(gr.in_degree(i) >= 2 && gr.out_degree(i) >= 2) continue;
		decompose_trivial_vertex(i);
		cnt++;
	}
	return cnt;
}

int predecompose::init_super_edges()
{
	mev.clear();
	edge_iterator it1, it2;
	PEEI pei;
	for(pei = gr.edges(), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		vector<int> v;
		int s = (*it1)->source();
		v.push_back(s);
		mev.insert(PEV(*it1, v));
	}
	return 0;
}

int predecompose::init_inner_weights()
{
	edge_iterator it1, it2;
	PEEI pei;
	for(pei = gr.edges(), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		edge_descriptor e = (*it1);
		double w = gr.get_edge_weight(e);
		edge_info ei = gr.get_edge_info(e);
		ei.weight = w;
		gr.set_edge_info(e, ei);
	}
	return 0;
}

int predecompose::init_vertex_map()
{
	v2v.clear();
	for(int i = 0; i < gr.num_vertices(); i++)
	{
		v2v.push_back(i);
	}
	return 0;
}

int predecompose::init_nonzeroset()
{
	nonzeroset.clear();
	for(int i = 1; i < gr.num_vertices() - 1; i++)
	{
		if(gr.degree(i) <= 0) continue;
		nonzeroset.insert(i);
	}
	return 0;
}

int predecompose::decompose_vertex_replace(int root, MPID &pe2w)
{
	// reassign weights
	MID md;
	for(MPID::iterator it = pe2w.begin(); it != pe2w.end(); it++)
	{
		int e1 = it->first.first;
		int e2 = it->first.second;
		double w = it->second;
		if(md.find(e1) == md.end()) md.insert(PID(e1, w));
		else md[e1] += w;
		if(md.find(e2) == md.end()) md.insert(PID(e2, w));
		else md[e2] += w;
	}

	for(MID::iterator it = md.begin(); it != md.end(); it++)
	{
		edge_descriptor e = i2e[it->first];
		double w = it->second;
		gr.set_edge_weight(e, w);
	}

	for(MPID::iterator it = pe2w.begin(); it != pe2w.end(); it++)
	{
		int e1 = it->first.first;
		int e2 = it->first.second;
		double w = it->second;
		int e = merge_adjacent_edges(e1, e2, w);
	}
	assert(gr.degree(root) == 0);
	nonzeroset.erase(root);
	return 0;
}

int predecompose::decompose_trivial_vertex(int x)
{
	balance_vertex(x);

	MPID pe2w;
	edge_iterator it1, it2;
	edge_iterator ot1, ot2;
	PEEI pe1, pe2;
	for(pe1 = gr.in_edges(x), it1 = pe1.first, it2 = pe1.second; it1 != it2; it1++)
	{
		int e1 = e2i[*it1];
		double w1 = gr.get_edge_weight(*it1);
		for(pe2 = gr.out_edges(x), ot1 = pe2.first, ot2 = pe2.second; ot1 != ot2; ot1++)
		{
			int e2 = e2i[*ot1];
			double w2 = gr.get_edge_weight(*ot1);
			double w = w1 <= w2 ? w1 : w2;

			pe2w.insert(PPID(PI(e1, e2), w));
		}
	}
	decompose_vertex_replace(x, pe2w);
	return 0;
}

int predecompose::merge_adjacent_equal_edges(int x, int y)
{
	if(i2e[x] == null_edge) return -1;
	if(i2e[y] == null_edge) return -1;

	edge_descriptor xx = i2e[x];
	edge_descriptor yy = i2e[y];

	int xs = (xx)->source();
	int xt = (xx)->target();
	int ys = (yy)->source();
	int yt = (yy)->target();

	if(xt != ys && yt != xs) return -1;
	if(yt == xs) return merge_adjacent_equal_edges(y, x);
	
	assert(xt == ys);

	edge_descriptor p = gr.add_edge(xs, yt);

	int n = i2e.size();
	i2e.push_back(p);
	assert(e2i.find(p) == e2i.end());
	e2i.insert(PEI(p, n));

	double wx0 = gr.get_edge_weight(xx);
	double wy0 = gr.get_edge_weight(yy);
	assert(fabs(wx0 - wy0) <= SMIN);

	int lx1 = gr.get_edge_info(xx).length;
	int ly1 = gr.get_edge_info(yy).length;
	int lxt = gr.get_vertex_info(xt).length;
	int lxy = lx1 + ly1 + lxt;

	gr.set_edge_weight(p, wx0 * 0.5 + wy0 * 0.5);
	gr.set_edge_info(p, edge_info(lxy));

	vector<int> v = mev[xx];
	v.insert(v.end(), mev[yy].begin(), mev[yy].end());

	if(mev.find(p) != mev.end()) mev[p] = v;
	else mev.insert(PEV(p, v));

	double sum1 = gr.get_in_weights(xt);
	double sum2 = gr.get_out_weights(xt);
	double sum = (sum1 + sum2) * 0.5;
	double r1 = gr.get_vertex_weight(xt) * (wx0 + wy0) * 0.5 / sum;
	double r2 = gr.get_vertex_weight(xt) - r1;
	gr.set_vertex_weight(xt, r2);

	assert(i2e[n] == p);
	assert(e2i.find(p) != e2i.end());
	assert(e2i[p] == n);
	assert(e2i[i2e[n]] == n);

	remove_edge(x);
	remove_edge(y);

	if(gr.in_degree(xt) == 0 && gr.out_degree(xt) == 0)
	{
		assert(gr.degree(xt) == 0);
		nonzeroset.erase(xt);
	}

	return n;
}

int predecompose::remove_edge(int e)
{
	edge_descriptor ee = i2e[e];
	assert(ee != null_edge);
	int s = ee->source();
	int t = ee->target();

	e2i.erase(ee);
	i2e[e] = null_edge;
	gr.remove_edge(ee);

	return 0;
}

int predecompose::merge_adjacent_edges(int x, int y, double ww)
{
	if(i2e[x] == null_edge) return -1;
	if(i2e[y] == null_edge) return -1;

	edge_descriptor xx = i2e[x];
	edge_descriptor yy = i2e[y];

	int xs = xx->source();
	int xt = xx->target();
	int ys = yy->source();
	int yt = yy->target();

	if(xt != ys) return merge_adjacent_edges(y, x, ww);
	assert(xt == ys);

	int x1 = split_edge(x, ww);
	int y1 = split_edge(y, ww);
	int xy = merge_adjacent_equal_edges(x1, y1);

	return xy;
}

int predecompose::merge_adjacent_edges(int x, int y)
{

	if(i2e[x] == null_edge) return -1;
	if(i2e[y] == null_edge) return -1;

	edge_descriptor xx = i2e[x];
	edge_descriptor yy = i2e[y];

	double wx = gr.get_edge_weight(xx);
	double wy = gr.get_edge_weight(yy);
	double ww = (wx <= wy) ? wx : wy;

	return merge_adjacent_edges(x, y, ww);
}

int predecompose::split_edge(int ei, double w)
{
	assert(i2e[ei] != null_edge);
	edge_descriptor ee = i2e[ei];

	double ww = gr.get_edge_weight(ee);

	if(fabs(ww - w) <= SMIN) return ei;
	assert(ww >= w + SMIN);

	int s = ee->source();
	int t = ee->target();

	edge_descriptor p2 = gr.add_edge(s, t);
	edge_info eif = gr.get_edge_info(ee);

	gr.set_edge_weight(ee, ww - w);		// old edge
	gr.set_edge_info(ee, eif);			// old edge
	gr.set_edge_weight(p2, w);			// new edge
	gr.set_edge_info(p2, eif);			// new edge

	if(mev.find(p2) != mev.end()) mev[p2] = mev[ee];
	else mev.insert(PEV(p2, mev[ee]));

	int n = i2e.size();
	i2e.push_back(p2);
	e2i.insert(PEI(p2, n));

	return n;
}

int predecompose::balance_vertex(int v)
{
	if(gr.degree(v) <= 0) return 0;

	edge_iterator it1, it2;
	PEEI pei;
	double w1 = 0, w2 = 0;
	for(pei = gr.in_edges(v), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		double w = gr.get_edge_weight(*it1);
		w1 += w;
	}
	for(pei = gr.out_edges(v), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		double w = gr.get_edge_weight(*it1);
		w2 += w;
	}

	assert(w1 >= SMIN);
	assert(w2 >= SMIN);

	// use max-meature
	//double ww = (wv >= w1 && wv >= w2) ? wv : (w1 >= w2 ? w1 : w2);
	//assert(ww >= w1 && ww >= w2);

	// use sqrt-meature
	double ww = sqrt(w1 * w2);

	// use convex combination
	//double ww = sqrt(0.5 * w1 * w1 + 0.5 * w2 * w2);

	double r1 = ww / w1;
	double r2 = ww / w2;

	double m1 = 0, m2 = 0;
	for(pei = gr.in_edges(v), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		double wx = gr.get_edge_weight(*it1);
		double wy = wx * r1;
		if(wy < 1.0)
		{
			m1 += (1.0 - wy);
			wy = 1.0;
		}
		gr.set_edge_weight(*it1, wy);
	}
	for(pei = gr.out_edges(v), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		double wx = gr.get_edge_weight(*it1);
		double wy = wx * r2;
		if(wy < 1.0)
		{
			m2 += 1.0 - wy;
			wy = 1.0;
		}
		gr.set_edge_weight(*it1, wy);
	}

	if(m1 > m2)
	{
		edge_descriptor e = gr.max_out_edge(v);
		double w = gr.get_edge_weight(e);
		gr.set_edge_weight(e, w + m1 - m2);
	}
	else if(m1 < m2)
	{
		edge_descriptor e = gr.max_in_edge(v);
		double w = gr.get_edge_weight(e);
		gr.set_edge_weight(e, w + m2 - m1);
	}

	return 0;
}

int predecompose::collect_paths()
{
	for(int i = 0; i < i2e.size(); i++)
	{
		if(i2e[i] == null_edge) continue;
		int s = i2e[i]->source();
		int t = i2e[i]->target();
		collect_path(i, s, t);
	}
	return 0;
}

int predecompose::collect_path(int e, int s, int t)
{
	assert(mev.find(i2e[e]) != mev.end());

	vector<int> v0 = mev[i2e[e]];
	vector<int> v;
	for(int i = 0; i < v0.size(); i++) 
	{
		if(v2v[v0[i]] < 0) continue;
		v.push_back(v2v[v0[i]]);
	}

	sort(v.begin(), v.end());
	assert(v.front() == s);
	assert(v.back() < t);
	v.push_back(t);

	path p;
	p.abd = gr.get_edge_weight(i2e[e]);
	p.v = v;
	//p.v.clear();
	//p.v.assign(v.begin() + 1, v.end());
	paths.push_back(p);

	gr.remove_edge(i2e[e]);
	e2i.erase(i2e[e]);
	i2e[e] = null_edge;

	return 0;
}

int predecompose::write(ostream &os)
{
	os<<fixed;
	os.precision(2);
	for(int i = 0; i < paths.size(); i++)
	{
		const vector<int> &v = paths[i].v;
		int n = v.size();
		if(n <= 2) continue;

		double w = paths[i].abd;
		os << "pred " << n;
		for(int k = 0; k < n; k++) os << " " << v[k];
		os << " " << w << " 1";
		os << endl;
	}
	return 0;
}
