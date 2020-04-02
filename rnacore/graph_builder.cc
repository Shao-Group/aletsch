#include <cassert>
#include <cstdio>
#include <map>
#include <iomanip>
#include <fstream>

#include "constants.h"
#include "graph_builder.h"
#include "region.h"
#include "util.h"
#include "undirected_graph.h"

graph_builder::graph_builder(bundle &b)
	: bd(b)
{}

int graph_builder::build(splice_graph &gr)
{
	build_junctions();
	build_regions();
	build_partial_exons();
	link_partial_exons();
	build_splice_graph(gr);
	//revise_splice_graph_full(gr, cfg);

	return 0;
}

int graph_builder::clear()
{
	junctions.clear();
	regions.clear();
	pexons.clear();
	return 0;
}

int graph_builder::build_junctions()
{
	map<PI32, vector<int> > m;
	for(int i = 0; i < bd.hits.size(); i++)
	{
		vector<int32_t> v = bd.hits[i].spos;
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
		//if(v.size() < cfg->min_splice_boundary_bd.hits) continue; TODO

		int32_t p1 = it->first.first;
		int32_t p2 = it->first.second;

		int s0 = 0;
		int s1 = 0;
		int s2 = 0;
		int nm = 0;
		for(int k = 0; k < v.size(); k++)
		{
			const hit &h = bd.hits[v[k]];
			nm += h.nm;
			if(h.xs == '.') s0++;
			if(h.xs == '+') s1++;
			if(h.xs == '-') s2++;
		}

		//printf("junction: %s:%d-%d (%d, %d, %d) %d\n", chrm.c_str(), p1, p2, s0, s1, s2, s1 < s2 ? s1 : s2);

		junction jc(p1, p2, v.size());
		jc.nm = nm;
		if(s1 == 0 && s2 == 0) jc.strand = '.';
		else if(s1 >= 1 && s2 >= 1) jc.strand = '.';
		else if(s1 > s2) jc.strand = '+';
		else jc.strand = '-';
		junctions.push_back(jc);

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

int graph_builder::build_regions()
{
	MPI s;
	s.insert(PI(bd.lpos, START_BOUNDARY));
	s.insert(PI(bd.rpos, END_BOUNDARY));
	for(int i = 0; i < junctions.size(); i++)
	{
		junction &jc = junctions[i];

		double ave, dev;
		evaluate_rectangle(bd.mmap, jc.lpos, jc.rpos, ave, dev);

		int32_t l = jc.lpos;
		int32_t r = jc.rpos;

		if(s.find(l) == s.end()) s.insert(PI(l, LEFT_SPLICE));
		else if(s[l] == RIGHT_SPLICE) s[l] = LEFT_RIGHT_SPLICE;

		if(s.find(r) == s.end()) s.insert(PI(r, RIGHT_SPLICE));
		else if(s[r] == LEFT_SPLICE) s[r] = LEFT_RIGHT_SPLICE;
	}

	for(int i = 0; i < pexons.size(); i++)
	{
		partial_exon &p = pexons[i];
		if(s.find(p.lpos) != s.end()) s.insert(PI(p.lpos, p.ltype));
		if(s.find(p.rpos) != s.end()) s.insert(PI(p.rpos, p.rtype));
	}

	vector<PPI> v(s.begin(), s.end());
	sort(v.begin(), v.end());

	regions.clear();
	for(int k = 0; k < v.size() - 1; k++)
	{
		int32_t l = v[k].first;
		int32_t r = v[k + 1].first;
		int ltype = v[k].second; 
		int rtype = v[k + 1].second; 

		if(ltype == LEFT_RIGHT_SPLICE) ltype = RIGHT_SPLICE;
		if(rtype == LEFT_RIGHT_SPLICE) rtype = LEFT_SPLICE;

		// TODO parameter
		//regions.push_back(region(l, r, ltype, rtype, &mmap, &imap, cfg->min_subregion_gap, cfg->min_subregion_length, cfg->min_subregion_overlap));
		regions.push_back(region(l, r, ltype, rtype, &(bd.mmap), &(bd.imap), 3, 15, 1.5));
	}

	return 0;
}

int graph_builder::build_partial_exons()
{
	pexons.clear();
	for(int i = 0; i < regions.size(); i++)
	{
		region &r = regions[i];
		pexons.insert(pexons.end(), r.pexons.begin(), r.pexons.end());
	}
	return 0;
}

int graph_builder::link_partial_exons()
{
	if(pexons.size() == 0) return 0;

	MPI lm;
	MPI rm;
	for(int i = 0; i < pexons.size(); i++)
	{
		int32_t l = pexons[i].lpos;
		int32_t r = pexons[i].rpos;

		assert(lm.find(l) == lm.end());
		assert(rm.find(r) == rm.end());
		lm.insert(PPI(l, i));
		rm.insert(PPI(r, i));
	}

	for(int i = 0; i < junctions.size(); i++)
	{
		junction &b = junctions[i];
		MPI::iterator li = rm.find(b.lpos);
		MPI::iterator ri = lm.find(b.rpos);

		assert(li != rm.end());
		assert(ri != lm.end());

		if(li != rm.end() && ri != lm.end())
		{
			b.lexon = li->second;
			b.rexon = ri->second;
		}
		else
		{
			b.lexon = b.rexon = -1;
		}
	}
	return 0;
}

int graph_builder::build_splice_graph(splice_graph &gr)
{
	gr.clear();
	gr.strand = bd.strand;
	gr.chrm = bd.chrm;

	// vertices: start, each region, end
	gr.add_vertex();
	vertex_info vi0;
	vi0.lpos = bd.lpos;
	vi0.rpos = bd.lpos;
	gr.set_vertex_weight(0, 0);
	gr.set_vertex_info(0, vi0);
	for(int i = 0; i < pexons.size(); i++)
	{
		const partial_exon &r = pexons[i];
		int length = r.rpos - r.lpos;
		assert(length >= 1);
		gr.add_vertex();
		gr.set_vertex_weight(i + 1, r.ave < 1.0 ? 1.0 : r.ave);
		vertex_info vi;
		vi.lpos = r.lpos;
		vi.rpos = r.rpos;
		vi.length = length;
		vi.stddev = r.dev;// < 1.0 ? 1.0 : r.dev;
		gr.set_vertex_info(i + 1, vi);
	}

	gr.add_vertex();
	vertex_info vin;
	vin.lpos = bd.rpos;
	vin.rpos = bd.rpos;
	gr.set_vertex_weight(pexons.size() + 1, 0);
	gr.set_vertex_info(pexons.size() + 1, vin);

	// edges: each junction => and e2w
	for(int i = 0; i < junctions.size(); i++)
	{
		const junction &b = junctions[i];

		if(b.lexon < 0 || b.rexon < 0) continue;

		const partial_exon &x = pexons[b.lexon];
		const partial_exon &y = pexons[b.rexon];

		edge_descriptor p = gr.add_edge(b.lexon + 1, b.rexon + 1);
		assert(b.count >= 1);
		edge_info ei;
		ei.weight = b.count;
		ei.strand = b.strand;
		gr.set_edge_info(p, ei);
		gr.set_edge_weight(p, b.count);
	}

	// edges: connecting start/end and pexons
	int ss = 0;
	int tt = pexons.size() + 1;
	for(int i = 0; i < pexons.size(); i++)
	{
		const partial_exon &r = pexons[i];

		if(r.ltype == START_BOUNDARY)
		{
			edge_descriptor p = gr.add_edge(ss, i + 1);
			double w = r.ave;
			if(i >= 1 && pexons[i - 1].rpos == r.lpos) w -= pexons[i - 1].ave;
			if(w < 1.0) w = 1.0;
			gr.set_edge_weight(p, w);
			edge_info ei;
			ei.weight = w;
			gr.set_edge_info(p, ei);
		}

		if(r.rtype == END_BOUNDARY) 
		{
			edge_descriptor p = gr.add_edge(i + 1, tt);
			double w = r.ave;
			if(i < pexons.size() - 1 && pexons[i + 1].lpos == r.rpos) w -= pexons[i + 1].ave;
			if(w < 1.0) w = 1.0;
			gr.set_edge_weight(p, w);
			edge_info ei;
			ei.weight = w;
			gr.set_edge_info(p, ei);
		}
	}

	// edges: connecting adjacent pexons => e2w
	for(int i = 0; i < (int)(pexons.size()) - 1; i++)
	{
		const partial_exon &x = pexons[i];
		const partial_exon &y = pexons[i + 1];

		if(x.rpos != y.lpos) continue;

		assert(x.rpos == y.lpos);
		
		int xd = gr.out_degree(i + 1);
		int yd = gr.in_degree(i + 2);
		double wt = (xd < yd) ? x.ave : y.ave;
		//int32_t xr = compute_overlap(mmap, x.rpos - 1);
		//int32_t yl = compute_overlap(mmap, y.lpos);
		//double wt = xr < yl ? xr : yl;

		edge_descriptor p = gr.add_edge(i + 1, i + 2);
		double w = (wt < 1.0) ? 1.0 : wt;
		gr.set_edge_weight(p, w);
		edge_info ei;
		ei.weight = w;
		gr.set_edge_info(p, ei);
	}

	return 0;
}

int graph_builder::print(int index)
{
	printf("graph-builder%d: ", index);

	// print regions
	for(int i = 0; i < regions.size(); i++)
	{
		regions[i].print(i);
	}

	// print junctions 
	for(int i = 0; i < junctions.size(); i++)
	{
		junctions[i].print(bd.chrm, i);
	}

	// print partial exons
	for(int i = 0; i < pexons.size(); i++)
	{
		pexons[i].print(i);
	}

	printf("\n");

	return 0;
}
