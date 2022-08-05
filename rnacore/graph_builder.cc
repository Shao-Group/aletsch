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

#include "constants.h"
#include "graph_builder.h"
#include "graph_reviser.h"
#include "region.h"
#include "util.h"
#include "undirected_graph.h"

graph_builder::graph_builder(bundle_base &b, const parameters &c, const sample_profile &s)
	: bd(b), cfg(c), sp(s)
{}

int graph_builder::build(splice_graph &gr)
{
	build_junctions();
	remove_opposite_junctions();
	build_regions();
	build_partial_exons();
	classify_partial_exons();
	link_partial_exons();
	build_splice_graph(gr);
	refine_splice_graph(gr);
	return 0;
}

int graph_builder::clear()
{
	junctions.clear();
	regions.clear();
	pexons.clear();
	regional.clear();
	return 0;
}

int graph_builder::build_junctions()
{
	chain_set jcst;

	for(int i = 0; i < bd.hcst.chains.size(); i++)
	{
		for(int j = 0; j < bd.hcst.chains[i].size(); j++)
		{
			vector<int32_t> &v = bd.hcst.chains[i][j].first;
			AI3 &a = bd.hcst.chains[i][j].second;

			if(v.size() <= 0) continue;
			if(v.size() % 2 != 0) continue;

			for(int k = 0; k < v.size() / 2; k++)
			{
				vector<int32_t> z;
				z.push_back(v[k * 2 + 0]);
				z.push_back(v[k * 2 + 1]);
				jcst.add(z, a);
			}
		}
	}

	for(int i = 0; i < bd.fcst.chains.size(); i++)
	{
		for(int j = 0; j < bd.fcst.chains[i].size(); j++)
		{
			vector<int32_t> &v = bd.fcst.chains[i][j].first;
			AI3 &a = bd.fcst.chains[i][j].second;

			if(v.size() <= 0) continue;
			if(v.size() % 2 != 0) continue;

			for(int k = 0; k < v.size() / 2; k++)
			{
				vector<int32_t> z;
				z.push_back(v[k * 2 + 0]);
				z.push_back(v[k * 2 + 1]);
				jcst.add(z, a);
			}
		}
	}

	for(int i = 0; i < jcst.chains.size(); i++)
	{
		for(int j = 0; j < jcst.chains[i].size(); j++)
		{
			vector<int32_t> &v = jcst.chains[i][j].first;
			AI3 &a = jcst.chains[i][j].second;

			if(v.size() != 2) continue;
			if(v[0] >= v[1]) continue;

			int count = a[0] + a[1] + a[2];
			if(count < cfg.min_junction_support) continue;

			junction jc(v[0], v[1], count);
			jc.xs0 = a[0];
			jc.xs1 = a[1];
			jc.xs2 = a[2];

			if(jc.xs1 > jc.xs2) jc.strand = '+';
			else if(jc.xs1 < jc.xs2) jc.strand = '-';
			else jc.strand = '.';

			//if(cfg.verbose >= 2) jc.print("chr1", i);

			junctions.push_back(jc);
		}
	}

	return 0;
}


int graph_builder::remove_opposite_junctions()
{
	double threshold = 15;
	set<int> fb;
	for(int i = 0; i < junctions.size(); i++)
	{
		if(fb.find(i) != fb.end()) continue;
		for(int j = i + 1; j < junctions.size(); j++)
		{
			if(fb.find(j) != fb.end()) continue;
			junction &x = junctions[i];
			junction &y = junctions[j];
			//if(x.strand == '.') continue;
			//if(y.strand == '.') continue;
			if(x.strand == y.strand) continue;

			double threshold = cfg.normal_junction_threshold;
			int32_t z = (x.rpos - x.lpos) - (y.rpos - y.lpos);
			if(z == 0 || x.lpos == y.lpos || x.rpos == y.rpos) threshold = cfg.extend_junction_threshold;

			double d = fabs(x.lpos - y.lpos) + fabs(x.rpos - y.rpos);
			if(d > threshold) continue;

			if(x.count > y.count && x.nm * 1.0 / x.count < y.nm * 1.0 / y.count) fb.insert(j);
			if(x.count < y.count && x.nm * 1.0 / x.count > y.nm * 1.0 / y.count) fb.insert(i);
		}
	}

	vector<junction> v;
	for(int i = 0; i < junctions.size(); i++)
	{
		if(fb.find(i) != fb.end()) 
		{
			if(cfg.verbose >= 2)
			{
				printf("remove opposite junction: ");
				junctions[i].print(bd.chrm, i);
			}
		}
		else
		{
			v.push_back(junctions[i]);
		}
	}

	junctions = v;
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

		//double ave, dev;
		//evaluate_rectangle(bd.mmap, jc.lpos, jc.rpos, ave, dev);

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

		regions.push_back(region(l, r, ltype, rtype, &(bd.mmap), &(bd.imap), cfg, sp));
	}

	return 0;
}

int graph_builder::build_partial_exons()
{
	pexons.clear();
	regional.clear();
	for(int i = 0; i < regions.size(); i++)
	{
		region &r = regions[i];
		for(int k = 0; k < r.pexons.size(); k++)
		{
			partial_exon &pe = r.pexons[k];
			pexons.push_back(pe);
			if((pe.lpos != bd.lpos || pe.rpos != bd.rpos) && pe.ltype == START_BOUNDARY && pe.rtype == END_BOUNDARY) regional.push_back(true);
			else regional.push_back(false);
		}
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

		if(ri == lm.end() || li == rm.end())
		{
			// print for testing
			printf("TEST: ");
			b.print("A", 999);
			for(int k = 0; k < regions.size(); k++) regions[k].print(k);
			for(int k = 0; k < pexons.size(); k++) pexons[k].print(k);
			for(int k = 0; k < junctions.size(); k++) junctions[k].print("X", k);
			printf("\n");
		}

		//if(li != rm.end()) printf("cannot find b.lpos = %d\n", b.lpos);
		//if(ri != lm.end()) printf("cannot find b.rpos = %d\n", b.rpos);

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

		//printf("link junction %d-%d to %d-%d\n", b.lpos, b.rpos, b.lexon, b.rexon);
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
	vi0.type = 0;
	gr.set_vertex_weight(0, 0);
	gr.set_vertex_info(0, vi0);
	for(int i = 0; i < pexons.size(); i++)
	{
		const partial_exon &r = pexons[i];
		int length = r.rpos - r.lpos;
		assert(length >= 1);
		gr.add_vertex();
		double w = r.ave;
		if(w < cfg.min_guaranteed_edge_weight) w = cfg.min_guaranteed_edge_weight;
		gr.set_vertex_weight(i + 1, w);
		vertex_info vi;
		if(r.pvalue < 0.5) vi.type = 0;
		else vi.type = 1;
		vi.lpos = r.lpos;
		vi.rpos = r.rpos;
		vi.stddev = r.dev;
		vi.maxcov = r.max;
		vi.length = length;
		vi.regional = regional[i];
		gr.set_vertex_info(i + 1, vi);
	}

	gr.add_vertex();
	vertex_info vin;
	vin.lpos = bd.rpos;
	vin.rpos = bd.rpos;
	vin.type = 0;
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
		if(b.strand == '+') ei.strand = 1;
		if(b.strand == '-') ei.strand = 2;
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
			if(w < cfg.min_guaranteed_edge_weight) w = cfg.min_guaranteed_edge_weight;
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
			if(w < cfg.min_guaranteed_edge_weight) w = cfg.min_guaranteed_edge_weight;
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

		double wt = x.ave;
		if(xd < yd) wt = x.ave;
		else if(xd > yd) wt = y.ave;
		else if(x.ave < y.ave) wt = x.ave;
		else if(x.ave > y.ave) wt = y.ave;
		//double wt = (xd < yd) ? x.ave : y.ave;
		//int32_t xr = compute_overlap(mmap, x.rpos - 1);
		//int32_t yl = compute_overlap(mmap, y.lpos);
		//double wt = xr < yl ? xr : yl;

		edge_descriptor p = gr.add_edge(i + 1, i + 2);
		if(wt < cfg.min_guaranteed_edge_weight) wt = cfg.min_guaranteed_edge_weight;
		gr.set_edge_weight(p, wt);
		edge_info ei;
		ei.weight = wt;
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

int graph_builder::analyze_junctions()
{
	double threshold = 15;
	for(int i = 0; i < junctions.size(); i++)
	{
		for(int j = i + 1; j < junctions.size(); j++)
		{
			junction &x = junctions[i];
			junction &y = junctions[j];
			if(x.strand == y.strand) continue;
			double d = fabs(x.lpos - y.lpos) + fabs(x.rpos - y.rpos);
			int32_t z = (x.rpos - x.lpos) - (y.rpos - y.lpos);
			if(z != 0 && d > threshold) continue;
			if(z == 0 && d > 2 * threshold) continue;
			x.print(bd.chrm, i);
			y.print(bd.chrm, j);
			printf("\n");
		}
	}
	return 0;
}

int graph_builder::classify_partial_exons()
{
	map<PI32, int> mj;
	for(int i = 0; i < junctions.size(); i++)
	{
		junction &jc = junctions[i];
		PI32 p(jc.lpos, jc.rpos);
		assert(mj.find(p) == mj.end());
		mj.insert(make_pair(p, i));
	}

	for(int i = 0; i < pexons.size(); i++)
	{
		partial_exon &pe = pexons[i];
		bool b = false;
		if(pe.lpos == bd.lpos) b = true;
		if(pe.rpos == bd.rpos) b = true;
		if(pe.ltype == RIGHT_SPLICE) b = true;
		if(pe.rtype == LEFT_SPLICE) b = true;
		if(pe.ltype == LEFT_SPLICE && pe.rtype == RIGHT_SPLICE)
		{
			PI32 p(pe.lpos, pe.rpos);
			if(mj.find(p) == mj.end()) b = true;
			else if(junctions[mj[p]].count < pe.ave) b = true;
		}

		if(b == true) 
		{
			pe.pvalue = 0;
		}
		else 
		{
			pe.pvalue = 1;
			//pe.ave *= 0.3;
		}
	}
	return 0;
}
