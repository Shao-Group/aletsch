#include <cassert>
#include <cstdio>
#include <map>
#include <iomanip>
#include <fstream>
#include <algorithm>

#include "graph_revise.h"
#include "config.h"
#include "undirected_graph.h"

int revise_splice_graph_full(splice_graph &gr, config *cfg)
{
	refine_splice_graph(gr);

	while(true)
	{
		bool b = false;

		b = remove_trivial_vertices(gr);
		if(b == true) continue;

		b = extend_boundaries(gr);
		if(b == true) continue;

		b = remove_inner_boundaries(gr);
		if(b == true) continue;

		b = remove_small_exons(gr, cfg->min_exon_length);
		if(b == true) refine_splice_graph(gr);
		if(b == true) continue;

		b = remove_small_junctions(gr);
		if(b == true) refine_splice_graph(gr);
		if(b == true) continue;

		b = keep_surviving_edges(gr, cfg->min_surviving_edge_weight);
		if(b == true) refine_splice_graph(gr);
		if(b == true) continue;

		b = remove_intron_contamination(gr, cfg->max_intron_contamination_coverage);
		if(b == true) continue;

		break;
	}

	refine_splice_graph(gr);

	return 0;
}

int revise_splice_graph(splice_graph &gr, config *cfg)
{
	refine_splice_graph(gr);

	while(true)
	{
		bool b = false;

		b = keep_surviving_edges(gr, cfg->min_surviving_edge_weight);
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
		if(gr.in_degree(i) != 1) continue;
		if(gr.out_degree(i) != 1) continue;

		PEEI pei = gr.in_edges(i);
		edge_iterator it1 = pei.first, it2 = pei.second;
		edge_descriptor e1 = (*it1);

		pei = gr.out_edges(i);
		it1 = pei.first;
		it2 = pei.second;
		edge_descriptor e2 = (*it1);
		vertex_info vi = gr.get_vertex_info(i);
		int s = e1->source();
		int t = e2->target();

		if(s != 0 && t != n) continue;
		if(s != 0 && gr.out_degree(s) == 1) continue;
		if(t != n && gr.in_degree(t) == 1) continue;

		if(vi.stddev >= 0.01) continue;

		//if(verbose >= 2) printf("remove inner boundary: vertex = %d, weight = %.2lf, length = %d, pos = %d-%d\n", i, gr.get_vertex_weight(i), vi.length, vi.lpos, vi.rpos);

		gr.clear_vertex(i);
		flag = true;
	}
	return flag;
}

bool remove_intron_contamination(splice_graph &gr, double ratio)
{
	bool flag = false;
	for(int i = 1; i < gr.num_vertices(); i++)
	{
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
		vertex_info vi = gr.get_vertex_info(i);

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
