/*
Part of aletsch
(c) 2020 by Mingfu Shao, The Pennsylvania State University
See LICENSE for licensing.
*/


#include <cassert>
#include <cstdio>
#include <map>
#include <iomanip>
#include <fstream>
#include <algorithm>

#include "graph_reviser.h"
#include "undirected_graph.h"
#include "parameters.h"
#include "essential.h"

int revise_splice_graph_full(splice_graph &gr, const parameters &cfg)
{
	refine_splice_graph(gr);

	while(true)
	{
		bool b = false;

		//b = remove_trivial_vertices(gr);
		//if(b == true) continue;

		b = extend_boundaries(gr);
		if(b == true) continue;

		b = remove_inner_boundaries(gr);
		if(b == true) continue;

		b = remove_small_exons(gr, cfg.min_exon_length);
		if(b == true) refine_splice_graph(gr);
		if(b == true) continue;

		b = remove_small_junctions(gr);
		if(b == true) refine_splice_graph(gr);
		if(b == true) continue;

		b = keep_surviving_edges(gr, cfg.min_surviving_edge_weight);
		if(b == true) refine_splice_graph(gr);
		if(b == true) continue;

		b = remove_intron_contamination(gr, cfg.max_intron_contamination_coverage);
		if(b == true) continue;

		break;
	}

	refine_splice_graph(gr);

	return 0;
}

int revise_splice_graph(splice_graph &gr, const parameters &cfg)
{
	refine_splice_graph(gr);

	while(true)
	{
		bool b = false;

		b = keep_surviving_edges(gr, cfg.min_surviving_edge_weight);
		if(b == true) refine_splice_graph(gr);
		if(b == true) continue;

		break;
	}

	refine_splice_graph(gr);

	return 0;
}

bool extend_boundaries(splice_graph &gr)
{
	edge_iterator it1, it2;
	PEEI pei;
	for(pei = gr.edges(), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		edge_descriptor e = (*it1);
		int s = e->source();
		int t = e->target();
		int32_t p = gr.get_vertex_info(t).lpos - gr.get_vertex_info(s).rpos;
		double we = gr.get_edge_weight(e);
		double ws = gr.get_vertex_weight(s);
		double wt = gr.get_vertex_weight(t);

		if(p <= 0) continue;
		if(s == 0) continue;
		if(t == gr.num_vertices() - 1) continue;

		bool b = false;
		if(gr.out_degree(s) == 1 && ws >= 10.0 * we * we + 10.0) b = true;
		if(gr.in_degree(t) == 1 && wt >= 10.0 * we * we + 10.0) b = true;

		if(b == false) continue;

		if(gr.out_degree(s) == 1)
		{
			edge_descriptor ee = gr.add_edge(s, gr.num_vertices() - 1);
			gr.set_edge_weight(ee, ws);
			gr.set_edge_info(ee, edge_info());
		}
		if(gr.in_degree(t) == 1)
		{
			edge_descriptor ee = gr.add_edge(0, t);
			gr.set_edge_weight(ee, wt);
			gr.set_edge_info(ee, edge_info());
		}

		gr.remove_edge(e);

		return true;
	}

	return false;
}

VE compute_maximal_edges(splice_graph &gr)
{
	typedef pair<double, edge_descriptor> PDE;
	vector<PDE> ve;

	undirected_graph ug;
	edge_iterator it1, it2;
	PEEI pei;
	for(int i = 0; i < gr.num_vertices(); i++) ug.add_vertex();
	for(pei = gr.edges(), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		edge_descriptor e = (*it1);
		double w = gr.get_edge_weight(e);
		int s = e->source();
		int t = e->target();
		if(s == 0) continue;
		if(t == gr.num_vertices() - 1) continue;
		//printf(" add edge to ug: s = %d, t = %d, gr.vs = %lu, ug.vs = %lu, w = %.2lf\n", s, t, gr.num_vertices(), ug.num_vertices(), w);
		ug.add_edge(s, t);
		ve.push_back(PDE(w, e));
	}

	vector<int> vv = ug.assign_connected_components();

	sort(ve.begin(), ve.end());

	for(int i = 1; i < ve.size(); i++) assert(ve[i - 1].first <= ve[i].first);

	VE x;
	set<int> sc;
	for(int i = ve.size() - 1; i >= 0; i--)
	{
		edge_descriptor e = ve[i].second;
		double w = gr.get_edge_weight(e);
		if(w < 1.5) break;
		int s = e->source();
		int t = e->target();
		if(s == 0) continue;
		if(t == gr.num_vertices()) continue;
		int c1 = vv[s];
		int c2 = vv[t];
		assert(c1 == c2);
		if(sc.find(c1) != sc.end()) continue;
		x.push_back(e);
		sc.insert(c1);
	}
	return x;
}

bool remove_trivial_vertices(splice_graph &gr)
{
	bool flag = false;
	for(int i = 1; i < gr.num_vertices() - 1; i++)
	{
		if(gr.in_degree(i) != 1) continue;
		if(gr.out_degree(i) != 1) continue;
		PEB p1 = gr.edge(0, i);
		PEB p2 = gr.edge(i, gr.num_vertices() - 1);
		if(p1.second == false) continue;
		if(p2.second == false) continue;
		gr.clear_vertex(i);
		flag = true;
	}
	return flag;
}

bool remove_small_exons(splice_graph &gr, int min_exon)
{
	bool flag = false;
	for(int i = 1; i < gr.num_vertices() - 1; i++)
	{
		vertex_info vi = gr.get_vertex_info(i);
		//if(vi.type == EMPTY_VERTEX) continue;

		bool b = true;
		edge_iterator it1, it2;
		PEEI pei;
		int32_t p1 = gr.get_vertex_info(i).lpos;
		int32_t p2 = gr.get_vertex_info(i).rpos;

		if(p2 - p1 >= min_exon) continue;
		if(gr.degree(i) <= 0) continue;

		for(pei = gr.in_edges(i), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
		{
			edge_descriptor e = (*it1);
			int s = e->source();
			//if(gr.out_degree(s) <= 1) b = false;
			if(s != 0 && gr.get_vertex_info(s).rpos == p1) b = false;
			if(b == false) break;
		}
		for(pei = gr.out_edges(i), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
		{
			edge_descriptor e = (*it1);
			int t = e->target();
			//if(gr.in_degree(t) <= 1) b = false;
			if(t != gr.num_vertices() - 1 && gr.get_vertex_info(t).lpos == p2) b = false;
			if(b == false) break;
		}

		if(b == false) continue;

		// only consider boundary small exons
		if(gr.edge(0, i).second == false && gr.edge(i, gr.num_vertices() - 1).second == false) continue;

		gr.clear_vertex(i);
		//vi.type = EMPTY_VERTEX;
		//gr.set_vertex_info(i, vi);

		flag = true;
	}
	return flag;
}

bool remove_small_junctions(splice_graph &gr)
{
	SE se;
	for(int i = 1; i < gr.num_vertices() - 1; i++)
	{
		if(gr.degree(i) <= 0) continue;

		bool b = true;
		edge_iterator it1, it2;
		PEEI pei;
		int32_t p1 = gr.get_vertex_info(i).lpos;
		int32_t p2 = gr.get_vertex_info(i).rpos;
		double wi = gr.get_vertex_weight(i);

		// compute max in-adjacent edge
		double ws = 0;
		for(pei = gr.in_edges(i), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
		{
			edge_descriptor e = (*it1);
			int s = e->source();
			double w = gr.get_vertex_weight(s);
			if(s == 0) continue;
			if(gr.get_vertex_info(s).rpos != p1) continue;
			if(w < ws) continue;
			ws = w;
		}

		// remove small in-junction
		for(pei = gr.in_edges(i), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
		{
			edge_descriptor e = (*it1);
			int s = e->source();
			double w = gr.get_edge_weight(e);
			if(s == 0) continue;
			if(gr.get_vertex_info(s).rpos == p1) continue;
			if(ws < 2.0 * w * w + 18.0) continue;
			if(wi < 2.0 * w * w + 18.0) continue;

			se.insert(e);
		}

		// compute max out-adjacent edge
		double wt = 0;
		for(pei = gr.out_edges(i), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
		{
			edge_descriptor e = (*it1);
			int t = e->target();
			double w = gr.get_vertex_weight(t);
			if(t == gr.num_vertices() - 1) continue;
			if(gr.get_vertex_info(t).lpos != p2) continue;
			if(w < wt) continue;
			wt = w;
		}

		// remove small in-junction
		for(pei = gr.out_edges(i), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
		{
			edge_descriptor e = (*it1);
			double w = gr.get_edge_weight(e);
			int t = e->target();
			if(t == gr.num_vertices() - 1) continue;
			if(gr.get_vertex_info(t).lpos == p2) continue;
			if(ws < 2.0 * w * w + 18.0) continue;
			if(wi < 2.0 * w * w + 18.0) continue;

			se.insert(e);
		}

	}

	if(se.size() <= 0) return false;

	for(SE::iterator it = se.begin(); it != se.end(); it++)
	{
		edge_descriptor e = (*it);
		gr.remove_edge(e);
	}

	return true;
}

bool remove_inner_boundaries(splice_graph &gr)
{
	bool flag = false;
	int n = gr.num_vertices() - 1;
	for(int i = 1; i < gr.num_vertices() - 1; i++)
	{
		vertex_info vi = gr.get_vertex_info(i);
		//if(vi.type == EMPTY_VERTEX) continue;

		if(gr.in_degree(i) != 1) continue;
		if(gr.out_degree(i) != 1) continue;

		PEEI pei = gr.in_edges(i);
		edge_iterator it1 = pei.first, it2 = pei.second;
		edge_descriptor e1 = (*it1);

		pei = gr.out_edges(i);
		it1 = pei.first;
		it2 = pei.second;
		edge_descriptor e2 = (*it1);
		int s = e1->source();
		int t = e2->target();

		if(s != 0 && t != n) continue;
		if(s != 0 && gr.out_degree(s) == 1) continue;
		if(t != n && gr.in_degree(t) == 1) continue;

		if(vi.stddev >= 0.01) continue;

		//if(verbose >= 2) printf("remove inner boundary: vertex = %d, weight = %.2lf, length = %d, pos = %d-%d\n", i, gr.get_vertex_weight(i), vi.length, vi.lpos, vi.rpos);

		gr.clear_vertex(i);
		//vi.type = EMPTY_VERTEX;
		//gr.set_vertex_info(i, vi);

		flag = true;
	}
	return flag;
}

bool remove_intron_contamination(splice_graph &gr, double ratio)
{
	bool flag = false;
	for(int i = 1; i < gr.num_vertices(); i++)
	{
		vertex_info vi = gr.get_vertex_info(i);
		//if(vi.type == EMPTY_VERTEX) continue;

		if(gr.in_degree(i) != 1) continue;
		if(gr.out_degree(i) != 1) continue;

		edge_iterator it1, it2;
		PEEI pei = gr.in_edges(i);
		it1 = pei.first;
		edge_descriptor e1 = (*it1);
		pei = gr.out_edges(i);
		it1 = pei.first;
		edge_descriptor e2 = (*it1);
		int s = e1->source();
		int t = e2->target();
		double wv = gr.get_vertex_weight(i);

		if(s == 0) continue;
		if(t == gr.num_vertices() - 1) continue;
		if(gr.get_vertex_info(s).rpos != vi.lpos) continue;
		if(gr.get_vertex_info(t).lpos != vi.rpos) continue;

		PEB p = gr.edge(s, t);
		if(p.second == false) continue;

		edge_descriptor ee = p.first;
		double we = gr.get_edge_weight(ee);

		if(wv > we) continue;
		if(wv > ratio) continue;

		//if(verbose >= 2) printf("clear intron contamination %d, weight = %.2lf, length = %d, edge weight = %.2lf\n", i, wv, vi.length, we); 

		gr.clear_vertex(i);
		//vi.type = EMPTY_VERTEX;
		//gr.set_vertex_info(i, vi);


		flag = true;
	}
	return flag;
}

bool keep_surviving_edges(splice_graph &gr, double surviving)
{
	set<int> sv1;
	set<int> sv2;
	SE se;
	edge_iterator it1, it2;
	PEEI pei;
	for(pei = gr.edges(), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		int s = (*it1)->source();
		int t = (*it1)->target();
		double w = gr.get_edge_weight(*it1);
		int c = gr.get_edge_info(*it1).count;
		int32_t p1 = gr.get_vertex_info(s).rpos;
		int32_t p2 = gr.get_vertex_info(t).lpos;
		if(w < surviving) continue;
		se.insert(*it1);
		sv1.insert(t);
		sv2.insert(s);
	}

	VE me = compute_maximal_edges(gr);
	for(int i = 0; i < me.size(); i++)
	{
		edge_descriptor ee = me[i];
		se.insert(ee);
		sv1.insert(ee->target());
		sv2.insert(ee->source());
	}

	while(true)
	{
		bool b = false;
		for(SE::iterator it = se.begin(); it != se.end(); it++)
		{
			edge_descriptor e = (*it);
			int s = e->source(); 
			int t = e->target();
			if(sv1.find(s) == sv1.end() && s != 0)
			{
				edge_descriptor ee = gr.max_in_edge(s);
				assert(ee != null_edge);
				assert(se.find(ee) == se.end());
				se.insert(ee);
				sv1.insert(s);
				sv2.insert(ee->source());
				b = true;
			}
			if(sv2.find(t) == sv2.end() && t != gr.num_vertices() - 1)
			{
				edge_descriptor ee = gr.max_out_edge(t);
				assert(ee != null_edge);
				assert(se.find(ee) == se.end());
				se.insert(ee);
				sv1.insert(ee->target());
				sv2.insert(t);
				b = true;
			}
			if(b == true) break;
		}
		if(b == false) break;
	}

	VE ve;
	for(pei = gr.edges(), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		if(se.find(*it1) != se.end()) continue;
		ve.push_back(*it1);
	}

	for(int i = 0; i < ve.size(); i++)
	{
		//if(verbose >= 2) printf("remove edge (%d, %d), weight = %.2lf\n", ve[i]->source(), ve[i]->target(), gr.get_edge_weight(ve[i]));
		gr.remove_edge(ve[i]);
	}

	if(ve.size() >= 1) return true;
	else return false;
}

int keep_surviving_edges(splice_graph &gr, const set<PI32> &js, double surviving)
{
	set<int> sv1;
	set<int> sv2;
	SE se;
	edge_iterator it1, it2;
	PEEI pei;
	for(pei = gr.edges(), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		int s = (*it1)->source();
		int t = (*it1)->target();
		double w = gr.get_edge_weight(*it1);
		int c = gr.get_edge_info(*it1).count;
		int32_t p1 = gr.get_vertex_info(s).rpos;
		int32_t p2 = gr.get_vertex_info(t).lpos;
		PI32 p(p1, p2);
		if(w < surviving && js.find(p) == js.end()) continue;
		se.insert(*it1);
		sv1.insert(t);
		sv2.insert(s);
	}

	while(true)
	{
		bool b = false;
		for(SE::iterator it = se.begin(); it != se.end(); it++)
		{
			edge_descriptor e = (*it);
			int s = e->source(); 
			int t = e->target();
			if(sv1.find(s) == sv1.end() && s != 0)
			{
				edge_descriptor ee = gr.max_in_edge(s);
				assert(ee != null_edge);
				assert(se.find(ee) == se.end());
				se.insert(ee);
				sv1.insert(s);
				sv2.insert(ee->source());
				b = true;
			}
			if(sv2.find(t) == sv2.end() && t != gr.num_vertices() - 1)
			{
				edge_descriptor ee = gr.max_out_edge(t);
				assert(ee != null_edge);
				assert(se.find(ee) == se.end());
				se.insert(ee);
				sv1.insert(ee->target());
				sv2.insert(t);
				b = true;
			}
			if(b == true) break;
		}
		if(b == false) break;
	}

	VE ve;
	for(pei = gr.edges(), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		if(se.find(*it1) != se.end()) continue;
		ve.push_back(*it1);
	}

	for(int i = 0; i < ve.size(); i++)
	{
		//if(verbose >= 2) printf("non-surviving edge (%d, %d), weight = %.2lf\n", ve[i]->source(), ve[i]->target(), gr.get_edge_weight(ve[i]));
		gr.remove_edge(ve[i]);
	}

	if(ve.size() >= 1) return true;
	else return false;
}

int keep_surviving_edges(splice_graph &gr, const set<int32_t> &js, double surviving)
{
	set<int> sv1;
	set<int> sv2;
	SE se;
	edge_iterator it1, it2;
	PEEI pei;
	int n = gr.num_vertices() - 1;
	for(pei = gr.edges(), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		int s = (*it1)->source();
		int t = (*it1)->target();
		double w = gr.get_edge_weight(*it1);
		int32_t p1 = gr.get_vertex_info(s).rpos;
		int32_t p2 = gr.get_vertex_info(t).lpos;
		//if(p1 == p2 && w < surviving) continue;
		//if(s == 0 && w < surviving) continue;
		//if(t == n && w < surviving) continue;
		if(p1 < p2 && w < surviving && (js.find(p1) == js.end() || js.find(p2) == js.end())) continue;
		se.insert(*it1);
		sv1.insert(t);
		sv2.insert(s);
	}

	while(true)
	{
		bool b = false;
		for(SE::iterator it = se.begin(); it != se.end(); it++)
		{
			edge_descriptor e = (*it);
			int s = e->source(); 
			int t = e->target();
			if(sv1.find(s) == sv1.end() && s != 0)
			{
				edge_descriptor ee = gr.max_in_edge(s);
				assert(ee != null_edge);
				assert(se.find(ee) == se.end());
				se.insert(ee);
				sv1.insert(s);
				sv2.insert(ee->source());
				b = true;
			}
			if(sv2.find(t) == sv2.end() && t != gr.num_vertices() - 1)
			{
				edge_descriptor ee = gr.max_out_edge(t);
				assert(ee != null_edge);
				assert(se.find(ee) == se.end());
				se.insert(ee);
				sv1.insert(ee->target());
				sv2.insert(t);
				b = true;
			}
			if(b == true) break;
		}
		if(b == false) break;
	}

	VE ve;
	for(pei = gr.edges(), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		if(se.find(*it1) != se.end()) continue;
		ve.push_back(*it1);
	}

	for(int i = 0; i < ve.size(); i++)
	{
		//if(verbose >= 2) printf("non-surviving edge (%d, %d), weight = %.2lf\n", ve[i]->source(), ve[i]->target(), gr.get_edge_weight(ve[i]));
		gr.remove_edge(ve[i]);
	}

	if(ve.size() >= 1) return true;
	else return false;
}


int keep_surviving_edges(splice_graph &gr, const set<int32_t> &js, const set<int32_t> &aj, double surviving)
{
	set<int> sv1;
	set<int> sv2;
	SE se;
	edge_iterator it1, it2;
	PEEI pei;
	for(pei = gr.edges(), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		int s = (*it1)->source();
		int t = (*it1)->target();
		double w = gr.get_edge_weight(*it1);
		int32_t p1 = gr.get_vertex_info(s).rpos;
		int32_t p2 = gr.get_vertex_info(t).lpos;
		if(p1 == p2 && w < surviving && aj.find(p1) == aj.end()) continue;
		if(p1 <  p2 && w < surviving && (js.find(p1) == js.end() || js.find(p2) == js.end())) continue;
		se.insert(*it1);
		sv1.insert(t);
		sv2.insert(s);
	}

	while(true)
	{
		bool b = false;
		for(SE::iterator it = se.begin(); it != se.end(); it++)
		{
			edge_descriptor e = (*it);
			int s = e->source(); 
			int t = e->target();
			if(sv1.find(s) == sv1.end() && s != 0)
			{
				edge_descriptor ee = gr.max_in_edge(s);
				assert(ee != null_edge);
				assert(se.find(ee) == se.end());
				se.insert(ee);
				sv1.insert(s);
				sv2.insert(ee->source());
				b = true;
			}
			if(sv2.find(t) == sv2.end() && t != gr.num_vertices() - 1)
			{
				edge_descriptor ee = gr.max_out_edge(t);
				assert(ee != null_edge);
				assert(se.find(ee) == se.end());
				se.insert(ee);
				sv1.insert(ee->target());
				sv2.insert(t);
				b = true;
			}
			if(b == true) break;
		}
		if(b == false) break;
	}

	VE ve;
	for(pei = gr.edges(), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		if(se.find(*it1) != se.end()) continue;
		ve.push_back(*it1);
	}

	for(int i = 0; i < ve.size(); i++)
	{
		//if(verbose >= 2) printf("non-surviving edge (%d, %d), weight = %.2lf\n", ve[i]->source(), ve[i]->target(), gr.get_edge_weight(ve[i]));
		gr.remove_edge(ve[i]);
	}

	if(ve.size() >= 1) return true;
	else return false;
}

int keep_surviving_edges(splice_graph &gr, const set<int32_t> &js, const set<int32_t> &aj, const set<int32_t> &sb, const set<int32_t> &tb, double surviving)
{
	set<int> sv1;
	set<int> sv2;
	SE se;
	edge_iterator it1, it2;
	PEEI pei;
	int n = gr.num_vertices() - 1;
	for(pei = gr.edges(), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		int s = (*it1)->source();
		int t = (*it1)->target();
		double w = gr.get_edge_weight(*it1);
		int32_t p1 = gr.get_vertex_info(s).rpos;
		int32_t p2 = gr.get_vertex_info(t).lpos;

		bool b = false;
		if(w >= surviving) b = true;
		else if(s == 0 && sb.find(p2) != sb.end()) b = true;
		else if(t == n && tb.find(p1) != tb.end()) b = true;
		else if(p1 == p2 && aj.find(p1) != aj.end()) b = true;
		else if(p1 <  p2 && js.find(p1) != js.end() && js.find(p2) != js.end()) b = true;

		if(b == false) continue;

		se.insert(*it1);
		sv1.insert(t);
		sv2.insert(s);
	}

	while(true)
	{
		bool b = false;
		for(SE::iterator it = se.begin(); it != se.end(); it++)
		{
			edge_descriptor e = (*it);
			int s = e->source(); 
			int t = e->target();
			if(sv1.find(s) == sv1.end() && s != 0)
			{
				edge_descriptor ee = gr.max_in_edge(s);
				assert(ee != null_edge);
				assert(se.find(ee) == se.end());
				se.insert(ee);
				sv1.insert(s);
				sv2.insert(ee->source());
				b = true;
			}
			if(sv2.find(t) == sv2.end() && t != gr.num_vertices() - 1)
			{
				edge_descriptor ee = gr.max_out_edge(t);
				assert(ee != null_edge);
				assert(se.find(ee) == se.end());
				se.insert(ee);
				sv1.insert(ee->target());
				sv2.insert(t);
				b = true;
			}
			if(b == true) break;
		}
		if(b == false) break;
	}

	VE ve;
	for(pei = gr.edges(), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		if(se.find(*it1) != se.end()) continue;
		ve.push_back(*it1);
	}

	for(int i = 0; i < ve.size(); i++)
	{
		printf("non-surviving edge (%d, %d), weight = %.2lf\n", ve[i]->source(), ve[i]->target(), gr.get_edge_weight(ve[i]));
		gr.remove_edge(ve[i]);
	}

	return 0;
}

int filter_start_boundaries(splice_graph &gr, const set<int32_t> &js, double surviving)
{
	VE ve;
	PEEI pei = gr.out_edges(0);
	int32_t z = gr.get_vertex_info(0).lpos;
	assert(z == gr.get_vertex_info(0).rpos);
	for(edge_iterator it = pei.first; it != pei.second; it++)
	{
		edge_descriptor e = (*it);
		double w = gr.get_edge_weight(e);
		if(w >= surviving) continue;

		int s = e->source();
		int t = e->target();
		assert(s == 0);

		int32_t p = gr.get_vertex_info(t).lpos;

		if(p == z) continue;
		if(js.find(p) != js.end()) continue;

		ve.push_back(e);
	}

	for(int i = 0; i < ve.size(); i++)
	{
		//if(verbose >= 2) printf("non-surviving start boundary (%d, %d), weight = %.2lf\n", ve[i]->source(), ve[i]->target(), gr.get_edge_weight(ve[i]));
		gr.remove_edge(ve[i]);
	}

	refine_splice_graph(gr);
	return 0;
}

int filter_end_boundaries(splice_graph &gr, const set<int32_t> &js, double surviving)
{
	VE ve;
	int n = gr.num_vertices() - 1;
	PEEI pei = gr.in_edges(n);
	int32_t z = gr.get_vertex_info(n).lpos;
	assert(z == gr.get_vertex_info(n).rpos);
	for(edge_iterator it = pei.first; it != pei.second; it++)
	{
		edge_descriptor e = (*it);
		double w = gr.get_edge_weight(e);
		if(w >= surviving) continue;

		int s = e->source();
		int t = e->target();

		int32_t p = gr.get_vertex_info(s).rpos;

		if(p == z) continue;
		if(js.find(p) != js.end()) continue;

		ve.push_back(e);
	}

	for(int i = 0; i < ve.size(); i++)
	{
		//if(verbose >= 2) printf("non-surviving end boundary (%d, %d), weight = %.2lf\n", ve[i]->source(), ve[i]->target(), gr.get_edge_weight(ve[i]));
		gr.remove_edge(ve[i]);
	}

	refine_splice_graph(gr);
	return 0;
}

int filter_junctions(splice_graph &gr, const set<int32_t> &js, double surviving)
{
	edge_iterator it1, it2;
	PEEI pei;
	VE ve;
	for(pei = gr.edges(), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		int s = (*it1)->source();
		int t = (*it1)->target();
		double w = gr.get_edge_weight(*it1);
		int32_t p1 = gr.get_vertex_info(s).rpos;
		int32_t p2 = gr.get_vertex_info(t).lpos;
		if(p1 <= p2) continue;
		if(w >= surviving) continue;
		if(js.find(p1) != js.end() && js.find(p2) == js.end()) continue;
		ve.push_back(*it1);
	}

	for(int i = 0; i < ve.size(); i++)
	{
		//if(verbose >= 2) printf("non-surviving junction (%d, %d), weight = %.2lf\n", ve[i]->source(), ve[i]->target(), gr.get_edge_weight(ve[i]));
		gr.remove_edge(ve[i]);
	}

	refine_splice_graph(gr);
	return 0;
}

int filter_graph(splice_graph &gr, const set<int32_t> &js, const set<int32_t> &aj, const set<int32_t> &sb, const set<int32_t> &tb, double surviving)
{
	VE ve;
	PEEI pei;
	edge_iterator it1, it2;
	int n = gr.num_vertices() - 1;
	for(pei = gr.edges(), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		int s = (*it1)->source();
		int t = (*it1)->target();
		double w = gr.get_edge_weight(*it1);
		int32_t p1 = gr.get_vertex_info(s).rpos;
		int32_t p2 = gr.get_vertex_info(t).lpos;

		bool b = false;
		if(w >= surviving) b = true;
		else if(s == 0 && sb.find(p2) != sb.end()) b = true;
		else if(t == n && tb.find(p1) != tb.end()) b = true;
		else if(p1 == p2 && aj.find(p1) != aj.end()) b = true;
		else if(p1 <  p2 && js.find(p1) != js.end() && js.find(p2) != js.end()) b = true;

		if(b == true) continue;
		ve.push_back(*it1);
	}

	for(int i = 0; i < ve.size(); i++)
	{
		printf("filter edge (%d, %d), weight = %.2lf\n", ve[i]->source(), ve[i]->target(), gr.get_edge_weight(ve[i]));
		gr.remove_edge(ve[i]);
	}

	refine_splice_graph(gr);
	return 0;
}

int refine_splice_graph(splice_graph &gr)
{
	while(true)
	{
		bool b = false;
		for(int i = 1; i < gr.num_vertices() - 1; i++)
		{
			if(gr.degree(i) == 0) continue;
			if(gr.in_degree(i) >= 1 && gr.out_degree(i) >= 1) continue;
			gr.clear_vertex(i);
			b = true;
		}
		if(b == false) break;
	}
	return 0;
}

int group_start_boundaries(splice_graph &gr, map<int32_t, int32_t> &smap, int32_t max_group_boundary_distance)
{
	smap.clear();
	vector<int> v;
	edge_iterator it;
	PEEI pei = gr.out_edges(0);
	for(edge_iterator it = pei.first; it != pei.second; it++)
	{
		edge_descriptor e = *it;
		assert(e->source() == 0);
		v.push_back(e->target());
	}

	if(v.size() <= 1) return 0;

	sort(v.begin(), v.end());

	int32_t p1 = gr.get_vertex_info(v[0]).lpos;
	int32_t p2 = p1;
	int k1 = v[0];
	int k2 = k1;
	PEB pa = gr.edge(0, v[0]);
	assert(pa.second == true);
	double wa = gr.get_edge_weight(pa.first);
	edge_info ea = gr.get_edge_info(pa.first);

	for(int i = 1; i < v.size(); i++)
	{
		int32_t p = gr.get_vertex_info(v[i]).lpos;
		PEB pb = gr.edge(0, v[i]);
		assert(pb.second == true);
		double wb = gr.get_edge_weight(pb.first);
		edge_info eb = gr.get_edge_info(pb.first);

		bool b = check_continuous_vertices(gr, k2, v[i]);

		assert(p >= p2);
		if(p - p2 > max_group_boundary_distance) b = false;

		if(b == false)
		{
			p1 = p;
			p2 = p;
			k1 = v[i];
			k2 = v[i];
			pa = pb;
			wa = wb;
			ea = eb;
		}
		else
		{
			smap.insert(pair<int32_t, int32_t>(p, p1));
			for(int j = k1; j < v[i]; j++)
			{
				PEB pc = gr.edge(j, j + 1);
				assert(pc.second == true);
				double vc = gr.get_vertex_weight(j);
				double wc = gr.get_edge_weight(pc.first);
				gr.set_vertex_weight(j, vc + wb);
				edge_info ec = gr.get_edge_info(pc.first);
				ec.count += eb.count;
				ec.weight += eb.weight;
				gr.set_edge_weight(pc.first, wc + wb);
				gr.set_edge_info(pc.first, ec);
			}
			wa += wb;
			ea.count += eb.count;
			ea.weight += eb.weight;
			gr.set_edge_weight(pa.first, wa);
			gr.set_edge_info(pa.first, ea);
			gr.remove_edge(pb.first);

			k2 = v[i];
			p2 = p;

			//if(verbose >= 2) printf("map start boundary %d:%d (weight = %.2lf) to %d:%d (weight = %.2lf)\n", v[i], p, wb, k1, p1, wa);
		}
	}
	return 0;
}

int group_end_boundaries(splice_graph &gr, map<int32_t, int32_t> &tmap, int32_t max_group_boundary_distance)
{
	tmap.clear();
	vector<int> v;
	edge_iterator it;
	PEEI pei = gr.in_edges(gr.num_vertices() - 1);
	for(edge_iterator it = pei.first; it != pei.second; it++)
	{
		edge_descriptor e = *it;
		assert(e->target() == gr.num_vertices() - 1);
		v.push_back(e->source());
	}

	if(v.size() <= 1) return 0;

	sort(v.begin(), v.end(), greater<int>());

	int32_t p1 = gr.get_vertex_info(v[0]).rpos;
	int32_t p2 = p1;
	int k1 = v[0];
	int k2 = k1;
	PEB pa = gr.edge(v[0], gr.num_vertices() - 1);
	assert(pa.second == true);
	double wa = gr.get_edge_weight(pa.first);

	for(int i = 1; i < v.size(); i++)
	{
		int32_t p = gr.get_vertex_info(v[i]).rpos;
		PEB pb = gr.edge(v[i], gr.num_vertices() - 1);
		assert(pb.second == true);
		double wb = gr.get_edge_weight(pb.first);

		bool b = check_continuous_vertices(gr, v[i], k2);

		assert(p <= p2);
		if(p2 - p > max_group_boundary_distance) b = false;

		if(b == false)
		{
			p1 = p;
			p2 = p;
			k1 = v[i];
			k2 = v[i];
			pa = pb;
			wa = wb;
		}
		else
		{
			tmap.insert(pair<int32_t, int32_t>(p, p1));
			for(int j = v[i]; j < k1; j++)
			{
				PEB pc = gr.edge(j, j + 1);
				assert(pc.second == true);
				double wc = gr.get_edge_weight(pc.first);
				double vc = gr.get_vertex_weight(j + 1);
				gr.set_edge_weight(pc.first, wc + wb);
				gr.set_vertex_weight(j + 1, wc + wb);
			}
			wa += wb;
			gr.set_edge_weight(pa.first, wa);
			gr.remove_edge(pb.first);

			k2 = v[i];
			p2 = p;

			//if(verbose >= 2) printf("map end boundary %d:%d (weight = %.2lf) to %d:%d (weight = %.2lf)\n", v[i], p, wb, k1, p1, wa);
		}
	}
	return 0;
}

int identify_boundaries(splice_graph &gr, const parameters &cfg)
{
	while(true)
	{
		bool b1 = identify_start_boundary(gr, cfg.min_boundary_log_ratio);
		bool b2 = identify_end_boundary(gr, cfg.min_boundary_log_ratio);
		if(b1 == false && b2 == false) break;
	}
	return 0;
}

bool identify_start_boundary(splice_graph &gr, double min_ratio)
{
	int besta = -1;
	int bestb = -1;
	double bestr = 0;
	double bestw = 0;
	for(int x = 1; x < gr.num_vertices() - 1; x++)
	{
		int a = left_continuous_extend(gr, x);
		int b = x;
		//int b = right_continuous_extend(gr, x);
		if(a < 0 || b < 0 || a > b) continue;

		double maxcov, in_sum;
		int f = determine_start_boundary(gr, a, b, maxcov, in_sum);
		if(f < 0) continue;

		double r = log(2 + maxcov) / log(2 + in_sum);
		if(r < bestr) continue;

		bestr = r;
		besta = a;
		bestb = b;
		bestw = maxcov - in_sum;
	}

	if(besta < 0) return false;
	if(bestr < min_ratio) return false;

	edge_descriptor ee = gr.add_edge(0, besta);
	gr.set_edge_weight(ee, bestw);
	gr.set_edge_info(ee, edge_info());

	//if(cfg.verbose >= 2) printf("identify start boundary at %d, weight = %.1lf, ratio = %.2lf, region = (%d, %d)\n", gr.get_vertex_info(besta).lpos, bestw, bestr, besta, bestb);

	return true;
}

bool identify_end_boundary(splice_graph &gr, double min_ratio)
{
	int besta = -1;
	int bestb = -1;
	double bestr = 0;
	double bestw = 0;
	for(int x = 1; x < gr.num_vertices() - 1; x++)
	{
		//int a = left_continuous_extend(gr, x);
		int a = x;
		int b = right_continuous_extend(gr, x);
		if(a < 0 || b < 0 || a > b) continue;

		double maxcov, out_sum;
		int f = determine_end_boundary(gr, a, b, maxcov, out_sum);
		if(f < 0) continue;

		double r = log(2 + maxcov) / log(2 + out_sum);
		if(r < bestr) continue;

		bestr = r;
		besta = a;
		bestb = b;
		bestw = maxcov - out_sum;
	}

	if(bestb < 0) return false;
	if(bestr < min_ratio) return false;

	edge_descriptor ee = gr.add_edge(bestb, gr.num_vertices() - 1);
	gr.set_edge_weight(ee, bestw);
	gr.set_edge_info(ee, edge_info());

	//printf("identify end boundary at %d, weight = %.1lf, ratio = %.2lf, region = (%d, %d)\n", gr.get_vertex_info(bestb).rpos, bestw, bestr, besta, bestb);

	return true;
}

int determine_start_boundary(splice_graph &gr, int a, int b, double &maxcov, double &in_sum)
{
	maxcov = 0;
	in_sum = 0;

	//gr.print(); printf("determine start boundary: a = %d, b = %d\n", a, b);
	for(int k = a; k <= b; k++)
	{
		if(gr.edge(0, k).second == true) return -1;
		//if(maxcov < gr.get_vertex_info(k).maxcov) maxcov = gr.get_vertex_info(k).maxcov;
		if(maxcov < gr.get_vertex_weight(k)) maxcov = gr.get_vertex_weight(k);

		PEEI pei = gr.in_edges(k);
		for(edge_iterator it = pei.first; it != pei.second; it++)
		{
			edge_descriptor e = *it;
			assert(e->target() == k);
			int v = e->source();
			if(v >= a && v <= b) assert(v == k - 1);
			if(v >= a && v <= b) continue;
			assert(v < a);
			assert(v != 0);
			in_sum += gr.get_edge_weight(e);
		}
	}
	return 0;
}

int determine_end_boundary(splice_graph &gr, int a, int b, double &maxcov, double &out_sum)
{
	maxcov = 0;
	out_sum = 0;
	int n = gr.num_vertices() - 1;
	for(int k = a; k <= b; k++)
	{
		if(gr.edge(k, n).second == true) return -1;
		//if(maxcov < gr.get_vertex_info(k).maxcov) maxcov = gr.get_vertex_info(k).maxcov;
		if(maxcov < gr.get_vertex_weight(k)) maxcov = gr.get_vertex_weight(k);

		PEEI pei = gr.out_edges(k);
		for(edge_iterator it = pei.first; it != pei.second; it++)
		{
			edge_descriptor e = *it;
			assert(e->source() == k);
			int v = e->target();
			if(v >= a && v <= b) assert(v == k + 1);
			if(v >= a && v <= b) continue;
			assert(v > b);
			assert(v != n);
			out_sum += gr.get_edge_weight(e);
		}
	}
	return 0;
}

int left_continuous_extend(splice_graph &gr, int x)
{
	int z = -1;				// will be returned
	int k = x;				// current vertex
	set<int> fb;
	while(true)
	{
		if(k <= 0) break;
		if(fb.find(k) != fb.end()) break;

		z = k;
		add_distant_in_vertices(gr, k, fb);

		if(k - 1 <= 0) break;
		if(gr.edge(k - 1, k).second == false) break;
		if(gr.get_vertex_info(k - 1).rpos != gr.get_vertex_info(k).lpos) break;

		k = k - 1;
	}
	return z;
}

int right_continuous_extend(splice_graph &gr, int x)
{
	int z = -1;				// will be returned
	int k = x;				// current vertex
	set<int> fb;
	int n = gr.num_vertices() - 1;
	while(true)
	{
		if(k >= n) break;
		if(fb.find(k) != fb.end()) break;

		z = k;
		add_distant_out_vertices(gr, k, fb);

		if(k + 1 >= n) break;
		if(gr.edge(k, k + 1).second == false) break;
		if(gr.get_vertex_info(k + 1).lpos != gr.get_vertex_info(k).rpos) break;

		k = k + 1;
	}
	return z;
}

int add_distant_in_vertices(splice_graph &gr, int x, set<int> &s)
{
	PEEI pei = gr.in_edges(x);
	int32_t p = gr.get_vertex_info(x).lpos;
	for(edge_iterator it = pei.first; it != pei.second; it++)
	{
		edge_descriptor e = *it;
		assert(e->target() == x);
		int v = e->source();
		if(gr.get_vertex_info(v).rpos == p) continue;
		s.insert(v);
	}
	return 0;
}

int add_distant_out_vertices(splice_graph &gr, int x, set<int> &s)
{
	PEEI pei = gr.out_edges(x);
	int32_t p = gr.get_vertex_info(x).rpos;
	for(edge_iterator it = pei.first; it != pei.second; it++)
	{
		edge_descriptor e = *it;
		assert(e->source() == x);
		int v = e->target();
		if(gr.get_vertex_info(v).lpos == p) continue;
		s.insert(v);
	}
	return 0;
}

int remove_false_boundaries(splice_graph &gr, bundle_base &bb, const parameters &cfg)
{
	map<int, int> fb1;		// end
	map<int, int> fb2;		// start
	for(int i = 0; i < bb.frgs.size(); i++)
	{
		if(bb.frgs[i][2] != 0) continue;
		hit &h1 = bb.hits[bb.frgs[i][0]];
		hit &h2 = bb.hits[bb.frgs[i][1]];
		int u1 = gr.locate_vertex(h1.rpos - 1);
		int u2 = gr.locate_vertex(h2.pos);

		if(u1 < 0 || u2 < 0) continue;
		if(u1 >= u2) continue;

		vertex_info v1 = gr.get_vertex_info(u1);
		vertex_info v2 = gr.get_vertex_info(u2);

		if(h1.rpos - v1.lpos <= cfg.bridge_end_relaxing) continue;
		if(v2.rpos - h2.pos  <= cfg.bridge_end_relaxing) continue;

		//printf("%s: u1 = %d, %d-%d, u2 = %d, %d-%d, h1.rpos = %d, h2.lpos = %d\n", h1.qname.c_str(), u1, v1.lpos, v1.rpos, u2, v2.lpos, v2.rpos, h1.rpos, h2.pos);

		if(gr.get_vertex_info(u1).rpos == h1.rpos)
		{
			if(fb1.find(u1) != fb1.end()) fb1[u1]++;
			else fb1.insert(make_pair(u1, 1));
		}

		if(gr.get_vertex_info(u2).lpos == h2.pos)
		{
			if(fb2.find(u2) != fb2.end()) fb2[u2]++;
			else fb2.insert(make_pair(u2, 1));
		}
	}

	for(auto &x : fb1)
	{
		PEB p = gr.edge(x.first, gr.num_vertices() - 1);
		vertex_info vi = gr.get_vertex_info(x.first);
		if(p.second == false) continue;
		double w = gr.get_vertex_weight(x.first);
		double z = log(1 + w) / log(1 + x.second);
		double s = log(1 + w) - log(1 + x.second);
		if(cfg.verbose >= 2) printf("detect false end boundary %d with %d reads, vertex = %d, w = %.2lf, type = %d, z = %.2lf, s = %.2lf\n", vi.rpos, x.second, x.first, w, vi.type, z, s); 
		if(s > 1) continue;
		vi.type = EMPTY_VERTEX;
		gr.set_vertex_info(x.first, vi);
		gr.remove_edge(p.first);
	}

	for(auto &x : fb2)
	{
		PEB p = gr.edge(0, x.first);
		vertex_info vi = gr.get_vertex_info(x.first);
		if(p.second == false) continue;
		double w = gr.get_vertex_weight(x.first);
		double z = log(1 + w) / log(1 + x.second);
		double s = log(1 + w) - log(1 + x.second);
		if(cfg.verbose >= 2) printf("detect false start boundary %d with %d reads, vertex = %d, w = %.2lf, type = %d, z = %.2lf, s = %.2lf\n", vi.lpos, x.second, x.first, w, vi.type, z, s); 
		if(s > 1) continue;
		vi.type = EMPTY_VERTEX;
		gr.set_vertex_info(x.first, vi);
		gr.remove_edge(p.first);
	}
	return 0;
}

int catch_false_boundaries(splice_graph &gr, bundle_base &bb, const parameters &cfg)
{
	map<PI, int> fb;		// end
	for(int i = 0; i < bb.frgs.size(); i++)
	{
		if(bb.frgs[i][2] != 0) continue;
		hit &h1 = bb.hits[bb.frgs[i][0]];
		hit &h2 = bb.hits[bb.frgs[i][1]];
		int u1 = gr.locate_vertex(h1.rpos - 1);
		int u2 = gr.locate_vertex(h2.pos);

		if(u1 < 0 || u2 < 0) continue;
		if(u1 >= u2) continue;

		vertex_info v1 = gr.get_vertex_info(u1);
		vertex_info v2 = gr.get_vertex_info(u2);
		if(h1.rpos - v1.lpos <= cfg.bridge_end_relaxing) continue;
		if(v2.rpos - h2.pos  <= cfg.bridge_end_relaxing) continue;


		//printf("%s: u1 = %d, %d-%d, u2 = %d, %d-%d, h1.rpos = %d, h2.lpos = %d\n", h1.qname.c_str(), u1, v1.lpos, v1.rpos, u2, v2.lpos, v2.rpos, h1.rpos, h2.pos);

		PI p(u1, u2);
		if(fb.find(p) != fb.end()) fb[p]++;
		else fb.insert(make_pair(p, 1));
	}

	for(auto &x : fb)
	{
		int u1 = x.first.first;
		int u2 = x.first.second;
		if(gr.check_path(u1, u2) == true) continue;
		double w1 = gr.get_vertex_weight(u1);
		double w2 = gr.get_vertex_weight(u2);
		double s1 = log(1 + w1) - log(1 + x.second);
		double s2 = log(1 + w2) - log(1 + x.second);
		if(cfg.verbose >= 2) printf("detect false end boundary: u1 = %d, u2 = %d, w1 = %.2lf, w2 = %.2lf, s1 = %.2lf, s2 = %.2lf, u1 = %d-%d, u2 = %d-%d\n",
				u1, u2, w1, w2, s1, s2, gr.get_vertex_info(u1).lpos, gr.get_vertex_info(u1).rpos, gr.get_vertex_info(u2).lpos, gr.get_vertex_info(u2).rpos); 
		//gr.remove_edge(p.first);
	}

	return 0;
}
