/*
Part of meta-scallop
(c) 2020 by Mingfu Shao, The Pennsylvania State University
See LICENSE for licensing.
*/

#include "bridge_solver.h"
#include "parameters.h"
#include "combined_graph.h"
#include "config.h"
#include "essential.h"
#include <sstream>
#include <algorithm>

combined_graph::combined_graph(const parameters &c)
	: cfg(c)
{
	num_combined = 0;
	strand = '?';
}

int combined_graph::copy_meta_information(const combined_graph &cb)
{
	sid = cb.sid;
	gid = cb.gid;
	chrm = cb.chrm;
	strand = cb.strand;
	return 0;
}

int combined_graph::set_gid(int batch, int instance, int subindex)
{
	char name[10240];
	sprintf(name, "instance.%d.%d.%d", batch, instance, subindex);
	gid = name;
	return 0;
}

int combined_graph::build(splice_graph &gr, const phase_set &p, const vector<pereads_cluster> &ub)
{
	chrm = gr.chrm;
	strand = gr.strand;
	num_combined = 1;

	build_regions(gr);
	build_start_bounds(gr);
	build_end_bounds(gr);
	build_splices_junctions(gr);
	ps = p;
	vc = ub;
	return 0;
}
	
int combined_graph::build_regions(splice_graph &gr)
{
	regions.clear();
	int n = gr.num_vertices() - 1;
	for(int i = 1; i < n; i++)
	{
		if(gr.degree(i) == 0) continue;
		double w = gr.get_vertex_weight(i);
		vertex_info vi = gr.get_vertex_info(i);
		PI32 p(vi.lpos, vi.rpos);
		regions.push_back(PPDI(p, DI(w, 1)));
	}
	return 0;
}

int combined_graph::build_start_bounds(splice_graph &gr)
{
	sbounds.clear();
	PEEI pei = gr.out_edges(0);
	int n = gr.num_vertices() - 1;
	for(edge_iterator it = pei.first; it != pei.second; it++)
	{
		int s = (*it)->source(); 
		int t = (*it)->target();
		assert(s == 0 && t > s);
		if(t == n) continue;
		int32_t p = gr.get_vertex_info(t).lpos;
		double w = gr.get_edge_weight(*it);
		PIDI pi(p, DI(w, 1));
		sbounds.push_back(pi);
	}
	return 0;
}

int combined_graph::build_end_bounds(splice_graph &gr)
{
	tbounds.clear();
	int n = gr.num_vertices() - 1;
	PEEI pei = gr.in_edges(n);
	for(edge_iterator it = pei.first; it != pei.second; it++)
	{
		int s = (*it)->source(); 
		int t = (*it)->target();
		assert(t == n);
		assert(s < t);
		if(s == 0) continue;
		int32_t p = gr.get_vertex_info(s).rpos;
		double w = gr.get_edge_weight(*it);
		PIDI pi(p, DI(w, 1));
		tbounds.push_back(pi);
	}
	return 0;
}

int combined_graph::build_splices_junctions(splice_graph &gr)
{
	junctions.clear();
	splices.clear();
	PEEI pei = gr.edges();
	set<int32_t> sp;
	int n = gr.num_vertices() - 1;
	for(edge_iterator it = pei.first; it != pei.second; it++)
	{
		int s = (*it)->source(); 
		int t = (*it)->target();
		assert(s < t);
		if(s == 0) continue;
		if(t == n) continue;
		double w = gr.get_edge_weight(*it);
		int32_t p1 = gr.get_vertex_info(s).rpos;
		int32_t p2 = gr.get_vertex_info(t).lpos;
		int c = 1;
		int strand = gr.get_edge_info(*it).strand;
		if(p1 >= p2) continue;

		TI32 p(PI32(p1, p2), strand);
		junctions.push_back(PTDI(p, DI(w, 1)));
		sp.insert(p1);
		sp.insert(p2);
	}
	splices.assign(sp.begin(), sp.end());
	sort(splices.begin(), splices.end());
	return 0;
}

int combined_graph::get_overlapped_splice_positions(const vector<int32_t> &v) const
{
	vector<int32_t> vv(v.size(), 0);
	vector<int32_t>::iterator it = set_intersection(v.begin(), v.end(), splices.begin(), splices.end(), vv.begin());
	return it - vv.begin();
}

int combined_graph::combine(combined_graph *cb)
{
	vector<combined_graph*> v;
	v.push_back(cb);
	combine(v);
	return 0;
}

int combined_graph::combine(vector<combined_graph*> &gv)
{
	if(gv.size() == 0) return 0;

	/*
	chrm = gv[0]->chrm;
	strand = gv[0]->strand;
	sp = gv[0]->sp;
	num_combined = 0;
	*/

	split_interval_double_map imap;
	map<TI32, DI> mj;
	map<int32_t, DI> ms;
	map<int32_t, DI> mt;

	combine_regions(imap);
	combine_junctions(mj);
	combine_start_bounds(ms);
	combine_end_bounds(mt);

	for(int i = 0; i < gv.size(); i++)
	{
		combined_graph *gt = gv[i];
		gt->combine_regions(imap);
		gt->combine_junctions(mj);
		gt->combine_start_bounds(ms);
		gt->combine_end_bounds(mt);
		ps.combine(gt->ps);
		num_combined += gt->num_combined;
	}

	regions.clear();
	for(SIMD it = imap.begin(); it != imap.end(); it++)
	{
		int32_t l = lower(it->first);
		int32_t r = upper(it->first);
		regions.push_back(PPDI(PI32(l, r), DI(it->second, 1)));
	}

	junctions.assign(mj.begin(), mj.end());
	sbounds.assign(ms.begin(), ms.end());
	tbounds.assign(mt.begin(), mt.end());

	return 0;
}

int combined_graph::combine_regions(split_interval_double_map &imap) const
{
	for(int i = 0; i < regions.size(); i++)
	{
		PI32 p = regions[i].first;
		double w = regions[i].second.first;
		imap += make_pair(ROI(p.first, p.second), w);
	}
	return 0;
}

int combined_graph::combine_junctions(map<TI32, DI> &m) const
{
	for(int i = 0; i < junctions.size(); i++)
	{
		TI32 p = junctions[i].first;
		DI d = junctions[i].second;

		map<TI32, DI>::iterator x = m.find(p);

		if(x == m.end())
		{
			m.insert(pair<TI32, DI>(p, d));
		}
		else 
		{
			x->second.first += d.first;
			x->second.second += d.second;
		}
	}
	return 0;
}

int combined_graph::combine_start_bounds(map<int32_t, DI> &m) const
{
	for(int i = 0; i < sbounds.size(); i++)
	{
		int32_t p = sbounds[i].first;
		DI d = sbounds[i].second;

		map<int32_t, DI>::iterator x = m.find(p);

		if(x == m.end())
		{
			m.insert(pair<int32_t, DI>(p, d));
		}
		else 
		{
			x->second.first += d.first;
			x->second.second += d.second;
		}
	}
	return 0;
}

int combined_graph::combine_end_bounds(map<int32_t, DI> &m) const
{
	for(int i = 0; i < tbounds.size(); i++)
	{
		int32_t p = tbounds[i].first;
		DI d = tbounds[i].second;

		map<int32_t, DI>::iterator x = m.find(p);

		if(x == m.end())
		{
			m.insert(pair<int32_t, DI>(p, d));
		}
		else 
		{
			x->second.first += d.first;
			x->second.second += d.second;
		}
	}
	return 0;
}

int combined_graph::append(const pereads_cluster &pc, const bridge_path &bbp)
{
	assert(bbp.type >= 0);
	append_regions(pc, bbp);
	append_junctions(pc, bbp);
	add_phases_from_bridged_pereads_cluster(pc, bbp, ps);
	return 0;
}

int combined_graph::append_regions(const pereads_cluster &pc, const bridge_path &bbp)
{
	int32_t p1, p2;
	if(bbp.chain.size() == 0)
	{
		p1 = pc.extend[1];
		p2 = pc.extend[2];
		if(p1 >= p2) return 0;
		PPDI pi(PI32(p1, p2), DI(pc.count, 1));
		regions.push_back(pi);
		return 0;
	}

	// append first region
	p1 = pc.extend[1];
	p2 = bbp.chain[0];
	if(p1 < p2)
	{
		PPDI pi(PI32(p1, p2), DI(pc.count, 1));
		regions.push_back(pi);
	}
	else
	{
		PPDI pi(PI32(p2, p1), DI(0.1, 1));
		regions.push_back(pi);
	}

	// append middle regions
	for(int i = 0; i < bbp.chain.size() / 2 - 1; i++)
	{
		p1 = bbp.chain[i * 2 + 1];
		p2 = bbp.chain[i * 2 + 2];
		assert(p1 < p2);
		PPDI pi(PI32(p1, p2), DI(pc.count, 1));
		regions.push_back(pi);
	}

	// append last region
	p1 = bbp.chain.back();
	p2 = pc.extend[2];
	if(p1 < p2)
	{
		PPDI pi(PI32(p1, p2), DI(pc.count, 1));
		regions.push_back(pi);
	}
	else
	{
		PPDI pi(PI32(p2, p1), DI(0.1, 1));
		regions.push_back(pi);
	}
	return 0;
}

int combined_graph::append_junctions(const pereads_cluster &pc, const bridge_path &bbp)
{
	// append middle regions
	for(int i = 0; i < bbp.chain.size() / 2; i++)
	{
		int32_t p1 = bbp.chain[i * 2 + 0];
		int32_t p2 = bbp.chain[i * 2 + 1];
		assert(p1 < p2);
		TI32 ti(PI32(p1, p2), bbp.strand);
		PTDI pi(ti, DI(pc.count, 1));
		junctions.push_back(pi);
	}
	return 0;
}

int combined_graph::build_splice_graph(splice_graph &gr, const parameters &cfg)
{
	gr.clear();

	gr.gid = gid;
	gr.chrm = chrm;
	gr.strand = strand;

	// add vertices
	gr.add_vertex();	// 0
	PIDI sb = get_leftmost_bound();
	vertex_info v0;
	v0.lpos = sb.first;
	v0.rpos = sb.first;
	gr.set_vertex_weight(0, 0);
	gr.set_vertex_info(0, v0);

	for(int i = 0; i < regions.size(); i++) 
	{
		gr.add_vertex();
		vertex_info vi;
		vi.lpos = regions[i].first.first;
		vi.rpos = regions[i].first.second;
		vi.count = regions[i].second.second;
		double w = regions[i].second.first;
		vi.length = vi.rpos - vi.lpos;
		gr.set_vertex_weight(i + 1, w);
		gr.set_vertex_info(i + 1, vi);
	}

	gr.add_vertex();	// n
	PIDI tb = get_rightmost_bound();
	vertex_info vn;
	vn.lpos = tb.first;
	vn.rpos = tb.first;
	gr.set_vertex_info(regions.size() + 1, vn);
	gr.set_vertex_weight(regions.size() + 1, 0);

	// build vertex index
	gr.build_vertex_index();

	// add sbounds
	for(int i = 0; i < sbounds.size(); i++)
	{
		int32_t p = sbounds[i].first;
		double w = sbounds[i].second.first;
		int c = sbounds[i].second.second;

		assert(gr.lindex.find(p) != gr.lindex.end());
		int k = gr.lindex[p];
		edge_descriptor e = gr.add_edge(0, k);
		edge_info ei;
		ei.weight = w;
		ei.count = c;
		gr.set_edge_info(e, ei);
		gr.set_edge_weight(e, w);
	}

	// add tbounds
	for(int i = 0; i < tbounds.size(); i++)
	{
		int32_t p = tbounds[i].first;
		double w = tbounds[i].second.first;
		int c = tbounds[i].second.second;

		assert(gr.rindex.find(p) != gr.rindex.end());
		int k = gr.rindex[p];
		edge_descriptor e = gr.add_edge(k, gr.num_vertices() - 1);
		edge_info ei;
		ei.weight = w;
		ei.count = c;
		gr.set_edge_info(e, ei);
		gr.set_edge_weight(e, w);
	}

	// shrink junctions first
	shrink_junctions();

	// add junctions
	for(int i = 0; i < junctions.size(); i++)
	{
		PI32 p = junctions[i].first.first;
		int strand = junctions[i].first.second;
		double w = junctions[i].second.first;
		int c = junctions[i].second.second;

		// TODO, should be asserted
		if(gr.rindex.find(p.first) == gr.rindex.end()) continue;
		if(gr.lindex.find(p.second) == gr.lindex.end()) continue;
		int s = gr.rindex[p.first];
		int t = gr.lindex[p.second];
		edge_descriptor e = gr.add_edge(s, t);
		edge_info ei;
		ei.weight = w;
		ei.count = c;
		ei.strand = strand;
		gr.set_edge_info(e, ei);
		gr.set_edge_weight(e, w);
	}

	// connect adjacent regions
	for(int i = 1; i < regions.size(); i++)
	{
		int32_t p1 = regions[i - 1].first.second;
		int32_t p2 = regions[i - 0].first.first;

		assert(p1 <= p2);
		if(p1 < p2) continue;

		PPDI ss = regions[i - 1];
		PPDI tt = regions[i - 0];
		if(ss.first.second != tt.first.first) continue;

		// TODO
		/*
		double w1 = gr.get_out_weights(i + 0);
		double w2 = gr.get_in_weights(i + 1);
		double ws = ss.second.first - w1;
		double wt = tt.second.first - w2;
		double w = (ws + wt) * 0.5;
		*/

		int xd = gr.out_degree(i + 0);
		int yd = gr.in_degree(i + 1);
		double w = (xd < yd) ? ss.second.first : tt.second.first;
		if(w < cfg.min_guaranteed_edge_weight) w = cfg.min_guaranteed_edge_weight;
		int c = ss.second.second;
		if(ss.second.second > tt.second.second) c = tt.second.second;
		edge_descriptor e = gr.add_edge(i + 0, i + 1);
		edge_info ei;
		ei.weight = w;
		ei.count = c;
		gr.set_edge_info(e, ei);
		gr.set_edge_weight(e, w);
	}
	return 0;
}

int combined_graph::shrink_junctions()
{
	vector< map<PI, int> > mm(3);
	set<PI> ss;
	for(int i = 0; i < junctions.size(); i++)
	{
		PI32 p = junctions[i].first.first;
		int s = junctions[i].first.second;
		assert(mm[s].find(p) == mm[s].end());
		mm[s].insert(make_pair(p, i));
		ss.insert(p);
	}

	vector<PTDI> v;
	for(auto &x : ss)
	{
		int c = -1;
		PTDI jx(make_pair(TI32(x, 0), DI(0, 0)));
		for(int k = 0; k < 3; k++)
		{
			if(mm[k].find(x) == mm[k].end()) continue;
			PTDI jc = junctions[mm[k][x]];
			assert(jc.first.first == jx.first.first);
			jx.second.first += jc.second.first;
			jx.second.second += jc.second.second;
			if(c == -1 || jc.second.first > junctions[mm[c][x]].second.first) c = k;
		}
		assert(c != -1);
		jx.first.second = c;
		v.push_back(jx);
	}
	junctions = v;
	return 0;
}

set<int32_t> combined_graph::get_reliable_splices(int samples, double weight)
{
	map<int32_t, DI> m;
	for(int i = 0; i < junctions.size(); i++)
	{
		PI32 p = junctions[i].first.first;
		int32_t p1 = p.first;
		int32_t p2 = p.second;
		double w = junctions[i].second.first;
		int c = junctions[i].second.second;

		if(m.find(p1) == m.end())
		{
			DI di(w, c);
			m.insert(pair<int32_t, DI>(p1, di));
		}
		else
		{
			m[p1].first += w;
			m[p1].second += c;
		}

		if(m.find(p2) == m.end())
		{
			DI di(w, c);
			m.insert(pair<int32_t, DI>(p2, di));
		}
		else
		{
			m[p2].first += w;
			m[p2].second += c;
		}
	}

	set<int32_t> s;
	for(map<int32_t, DI>::iterator it = m.begin(); it != m.end(); it++)
	{
		int32_t p = it->first;
		double w = it->second.first;
		int c = it->second.second;
		if(w < weight && c < samples) continue;
		s.insert(p);
	}
	return s;
}

int combined_graph::refine_junctions(vector<combined_graph*> &gv, const vector<sample_profile> &samples)
{
	map<TI32, int> mt;
	classify_junctions(gv, samples, mt);

	directed_graph gr;
	map<PI32, PI32> jm;
	build_junction_graph(gr, mt);
	build_junction_map(gr, jm);

	project_junctions(jm);
	for(int i = 0; i < gv.size(); i++) gv[i]->project_junctions(jm);

	return 0;
}

int combined_graph::classify_junctions(vector<combined_graph*> &gv, const vector<sample_profile> &samples, map<TI32, int> &mt)
{
	for(int i = 0; i < gv.size(); i++)
	{
		combined_graph *gt = gv[i];
		assert(gt->sid >= 0 && gt->sid < samples.size());

		for(int k = 0; k < gt->junctions.size(); k++)
		{
			PTDI p = gt->junctions[k];
			
			if(mt.find(p.first) == mt.end())
			{
				int type = 0;
				if(samples[gt->sid].data_type == "paired_end") type = 1;
				else if(samples[gt->sid].data_type == "ont") type = -1;
				else assert(false);
				mt.insert(make_pair(p.first, type));
			}
			else
			{
				mt[p.first] = 1;
			}
		}
	}
	return 0;
}

map<PI32, PI32> combined_graph::group_junctions()
{
	directed_graph gr;
	map<PI32, PI32> jm;
	build_junction_graph(gr);
	build_junction_map(gr, jm);
	return jm;
}

int combined_graph::project_junctions(const map<PI32, PI32> &jm)
{
	rebuild_junctions(jm);
	rebuild_splices();
	ps.project_junctions(jm);
	for(int i = 0; i < vc.size(); i++) vc[i].project_junctions(jm);
	return 0;
}

int combined_graph::rebuild_splices()
{
	set<int32_t> s;
	for(int i = 0; i < junctions.size(); i++)
	{
		s.insert(junctions[i].first.first.first);
		s.insert(junctions[i].first.first.second);
	}

	vector<int32_t> v(s.begin(), s.end());
	sort(v.begin(), v.end());
	splices = v;
	return 0;
}

int combined_graph::rebuild_junctions(const map<PI32, PI32> &jm)
{
	vector<PTDI> v;
	for(int i = 0; i < junctions.size(); i++)
	{
		PTDI z = junctions[i];
		PI32 p = z.first.first;
		map<PI32, PI32>::const_iterator it = jm.find(p);
		if(it != jm.end())
		{
			z.first.first.first = it->second.first;
			z.first.first.second = it->second.second;
		}
		v.push_back(z);
	}
	junctions = v;

	map<TI32, DI> mj;
	combine_junctions(mj);
	junctions.assign(mj.begin(), mj.end());

	return 0;
}

int combined_graph::build_junction_map(directed_graph &gr, map<PI32, PI32> &jm)
{
	set<int> fb;
	vector<int> topo = gr.topological_sort();
	for(int i = 0; i < topo.size(); i++)
	{
		int z = topo[i];
		if(fb.find(z) != fb.end()) continue;
		vector<int> v;
		vector<int> b;
		gr.bfs(z, v, b, fb);

		PI32 pz = junctions[z].first.first;

		// print
		if(v.size() >= 2)
		{
			printf("cluster junction weights: ");
			for(int k = 0; k < v.size(); k++)
			{
				int x = v[k];
				PTDI &p = junctions[x];
				printf("%.0lf, ", p.second.first);
			}
			printf("\n");
		}

		for(int k = 0; k < v.size(); k++)
		{
			int x = v[k];
			assert(fb.find(x) == fb.end());
			fb.insert(x);

			if(x == z) continue;

			PI32 px = junctions[x].first.first;
			if(jm.find(px) == jm.end()) jm.insert(make_pair(px, pz));

			//junctions[z].second.first += junctions[x].second.first;
			//junctions[x].second.first = -1;
		}

	}

	/*
	vector<PTDI> vv;
	for(int i = 0; i < junctions.size(); i++)
	{
		if(junctions[i].second.first < 0) continue;
		vv.push_back(junctions[i]);
	}
	junctions = vv;
	*/

	return 0;
}

int combined_graph::build_junction_graph(directed_graph &gr, const map<TI32, int> &mt)
{
	for(int i = 0; i < junctions.size(); i++)
	{
		gr.add_vertex();
	}

	for(int i = 0; i < junctions.size(); i++)
	{
		TI32 pi = junctions[i].first;
		map<TI32, int>::const_iterator xi = mt.find(pi);
		if(xi == mt.end()) continue;
		if(xi->second != 1) continue;
		for(int j = i + 1; j < junctions.size(); j++)
		{
			TI32 pj = junctions[j].first;
			map<TI32, int>::const_iterator xj = mt.find(pj);
			if(xj == mt.end()) continue;
			if(xi->second != -1) continue;
			int p = compare_two_junctions(junctions[i], junctions[j]);
			if(p == +1) gr.add_edge(i, j);
			//if(p == -1) gr.add_edge(j, i);
		}
	}
	return 0;
}

int combined_graph::build_junction_graph(directed_graph &gr)
{
	for(int i = 0; i < junctions.size(); i++)
	{
		gr.add_vertex();
	}

	for(int i = 0; i < junctions.size(); i++)
	{
		for(int j = i + 1; j < junctions.size(); j++)
		{
			int p = compare_two_junctions(junctions[i], junctions[j]);
			if(p == +1) gr.add_edge(i, j);
			if(p == -1) gr.add_edge(j, i);
		}
	}
	return 0;
}

int combined_graph::compare_two_junctions(PTDI &x, PTDI &y)
{
	//if(x.first.second != y.first.second) return 0;	// TODO
	PI &px = x.first.first;
	PI &py = y.first.first;

	int32_t dx = px.second - px.first;
	int32_t dy = py.second - py.first;
	if(fabs(dx - dy) > cfg.max_cluster_intron_distance) return 0;

	double s1 = fabs(px.first - py.first);
	double s2 = fabs(px.second - py.second);
	if(s1 > cfg.max_cluster_intron_shifting) return 0;
	if(s2 > cfg.max_cluster_intron_shifting) return 0;

	if(x.second.first > (y.second.first + 1) * (y.second.first + 1)) return +1;
	if(y.second.first > (x.second.first + 1) * (x.second.first + 1)) return -1;
	return 0;
}

int combined_graph::clear()
{
	num_combined = 0;
	sid = -1;
	gid = "";
	chrm = "";
	strand = '.';
	splices.clear();
	regions.clear();
	junctions.clear();
	sbounds.clear();
	tbounds.clear();
	ps.clear();
	vc.clear();
	return 0;
}

int combined_graph::print(int index)
{
	int pereads = 0;
	for(int i = 0; i < vc.size(); i++)
	{
		pereads += vc[i].chain1.size();
		pereads += vc[i].chain2.size();
		pereads += vc[i].bounds.size();
		pereads += vc[i].extend.size();
	}

	printf("combined-graph %d: sid = %d, gid = %s, #combined = %d, chrm = %s, strand = %c, #regions = %lu, #sbounds = %lu, #tbounds = %lu, #junctions = %lu, #phases = %lu, #pereads = %lu / %d\n", 
			index, sid, gid.c_str(), num_combined, chrm.c_str(), strand, regions.size(), sbounds.size(), tbounds.size(), junctions.size(), ps.pmap.size(), vc.size(), pereads);

	return 0;

	for(int i = 0; i < regions.size(); i++)
	{
		PI32 p = regions[i].first;
		DI d = regions[i].second;
		printf("region %d: [%d, %d), w = %.2lf, c = %d\n", i, p.first, p.second, d.first, d.second);
	}
	for(int i = 0; i < junctions.size(); i++)
	{
		TI32 p = junctions[i].first;
		DI d = junctions[i].second;
		printf("junction %d: [%d, %d, %d), w = %.2lf, c = %d\n", i, p.first.first, p.first.second, p.second, d.first, d.second);
	}
	for(int i = 0; i < sbounds.size(); i++)
	{
		int32_t p = sbounds[i].first;
		DI d = sbounds[i].second;
		printf("sbound %d: %d, w = %.2lf, c = %d\n", i, p, d.first, d.second);
	}
	for(int i = 0; i < tbounds.size(); i++)
	{
		int32_t p = tbounds[i].first;
		DI d = tbounds[i].second;
		printf("tbound %d: %d, w = %.2lf, c = %d\n", i, p, d.first, d.second);
	}
	ps.print();
	return 0;
}

PIDI combined_graph::get_leftmost_bound()
{
	PIDI x;
	x.first = -1;
	for(int i = 0; i < sbounds.size(); i++)
	{
		int32_t p = sbounds[i].first;
		if(x.first == -1 || p < x.first)
		{
			x = sbounds[i];
		}
	}
	return x;
}

PIDI combined_graph::get_rightmost_bound()
{
	PIDI x;
	x.first = -1;
	for(int i = 0; i < tbounds.size(); i++)
	{
		int32_t p = tbounds[i].first;
		if(x.first == -1 || p > x.first)
		{
			x = tbounds[i];
		}
	}
	return x;
}
