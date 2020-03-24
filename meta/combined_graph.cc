#include "meta_config.h"
#include "combined_graph.h"
#include "config.h"
#include "essential.h"
#include "graph_revise.h"
#include <sstream>
#include <algorithm>

combined_graph::combined_graph()
{
	num_combined = 0;
	strand = '?';
}

int combined_graph::build(splice_graph &gr, hyper_set &hs, vector<fcluster> &ub)
{
	chrm = gr.chrm;
	strand = gr.strand;
	num_combined = 1;

	build_regions(gr);
	build_start_bounds(gr);
	build_end_bounds(gr);
	build_splices_junctions(gr);
	build_phase(gr, hs);
	build_reads(gr, ub);
	return 0;
}
	
int combined_graph::build_regions(splice_graph &gr)
{
	regions.clear();
	int n = gr.num_vertices() - 1;
	for(int i = 1; i < n; i++)
	{
		double weight = gr.get_vertex_weight(i);
		vertex_info vi = gr.get_vertex_info(i);
		PI32 p(vi.lpos, vi.rpos);
		DI d(weight, 1);
		regions.push_back(PPDI(p, d));
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
		double w = gr.get_edge_weight(*it);
		int32_t p = gr.get_vertex_info(t).lpos;
		int c = 1;

		PIDI pi(p, DI(w, c));
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
		double w = gr.get_edge_weight(*it);
		int32_t p = gr.get_vertex_info(s).rpos;
		int c = 1;

		PIDI pi(p, DI(w, c));
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
		if(p1 >= p2) continue;

		PI32 p(p1, p2);
		DI d(w, c);
		junctions.push_back(PPDI(p, d));
		sp.insert(p1);
		sp.insert(p2);
	}
	splices.assign(sp.begin(), sp.end());
	sort(splices.begin(), splices.end());
	return 0;
}

int combined_graph::build_phase(splice_graph &gr, hyper_set &hs)
{
	phase.clear();
	map<vector<int32_t>, int> mm;
	for(MVII::const_iterator it = hs.nodes.begin(); it != hs.nodes.end(); it++)
	{
		const vector<int> &v = it->first;
		int w = it->second;
		int c = 1;

		if(v.size() <= 0) continue;
		vector<int32_t> vv;
		build_path_coordinates(gr, v, vv);

		if(vv.size() <= 1) continue;
		vector<int32_t> zz(vv.begin() + 1, vv.end() - 1);
		map<vector<int32_t>, int>::iterator tp = mm.find(zz);
		if(tp == mm.end()) 
		{
			vector<int32_t> pv;
			pv.push_back(vv.front());
			pv.push_back(vv.back());
			pv.push_back(w);
			mm.insert(pair<vector<int32_t>, int>(zz, phase.size()));
			phase.push_back(PV32(zz, pv));
		}
		else 
		{
			int k = tp->second;
			phase[k].second.push_back(vv.front());
			phase[k].second.push_back(vv.back());
			phase[k].second.push_back(w);
		}
	}
	return 0;
}

int combined_graph::build_reads(splice_graph &gr, vector<fcluster> &ub)
{
	reads.clear();
	map<vector<int32_t>, int> mm;
	int n = gr.num_vertices() - 1;
	for(int i = 0; i < ub.size(); i++)
	{
		fcluster &fc = ub[i];
		if(fc.v1.size() <= 0) continue;
		if(fc.v2.size() <= 0) continue;
		assert(fc.v1.front() != 0);
		assert(fc.v2.front() != 0);
		assert(fc.v1.back() != n);
		assert(fc.v2.back() != n);

		vector<int32_t> vv1;
		vector<int32_t> vv2;
		build_path_coordinates(gr, fc.v1, vv1);
		build_path_coordinates(gr, fc.v2, vv2);

		assert(vv1.size() >= 2);
		assert(vv2.size() >= 2);

		vector<int32_t> vv;
		vv.push_back(vv1.size() - 2);
		vv.push_back(vv2.size() - 2);
		vv.insert(vv.end(), vv1.begin() + 1, vv1.end() - 1);
		vv.insert(vv.end(), vv2.begin() + 1, vv2.end() - 1);

		vector<int32_t> uu;
		for(int j = 0; j < fc.frset.size(); j++)
		{
			fragment &fr = fc.frset[j];
			assert(vv1[1] - vv1[0] >= fr.k1l);
			assert(vv2[1] - vv2[0] >= fr.k2l);
			assert(vv1[vv1.size() - 1] - vv1[vv1.size() - 2] >= fr.k1r);
			assert(vv2[vv2.size() - 1] - vv2[vv2.size() - 2] >= fr.k2r);

			uu.push_back(vv1[0] + fr.k1l);
			uu.push_back(vv2[0] + fr.k2l);
			uu.push_back(vv1[vv1.size() - 1] - fr.k1r);
			uu.push_back(vv2[vv2.size() - 1] - fr.k2r);
		}

		map<vector<int32_t>, int>::iterator tp = mm.find(vv);
		if(tp == mm.end())
		{
			mm.insert(pair<vector<int32_t>, int>(vv, reads.size()));
			reads.push_back(PV32(vv, uu));
		}
		else
		{
			int k = tp->second;
			reads[k].second.insert(reads[k].second.end(), uu.begin(), uu.end());
		}
	}

	return 0;
}

int combined_graph::combine(const combined_graph &gt)
{
	if(children.size() == 0) children.push_back(*this);

	if(gt.children.size() == 0) children.push_back(gt);
	else children.insert(children.end(), gt.children.begin(), gt.children.end());

	if(chrm == "") chrm = gt.chrm;
	if(strand == '?') strand = gt.strand;
	assert(gt.chrm == chrm);
	assert(gt.strand == strand);

	num_combined += gt.num_combined;

	// combine splices
	vector<int32_t> vv(gt.splices.size() + splices.size(), 0);
	vector<int32_t>::iterator it = set_union(gt.splices.begin(), gt.splices.end(), splices.begin(), splices.end(), vv.begin());
	vv.resize(it - vv.begin());
	splices = vv;

	return 0;
}

int combined_graph::get_overlapped_splice_positions(const vector<int32_t> &v) const
{
	vector<int32_t> vv(v.size(), 0);
	vector<int32_t>::iterator it = set_intersection(v.begin(), v.end(), splices.begin(), splices.end(), vv.begin());
	return it - vv.begin();
}

int combined_graph::combine_children()
{
	if(children.size() == 0) return 0;

	split_interval_map imap;
	map<PI32, DI> mj;
	map<int32_t, DI> ms;
	map<int32_t, DI> mt;
	MV32 mp;
	MV32 mr;

	int num = 0;
	for(int i = 0; i < children.size(); i++)
	{
		combined_graph &gt = children[i];
		combine_regions(imap, gt);
		combine_junctions(mj, gt);
		combine_phase(mp, gt);
		combine_reads(mr, gt);
		combine_start_bounds(ms, gt);
		combine_end_bounds(mt, gt);
		num += gt.num_combined;
	}
	assert(num == num_combined);

	regions.clear();
	for(SIMI it = imap.begin(); it != imap.end(); it++)
	{
		int32_t l = lower(it->first);
		int32_t r = upper(it->first);
		regions.push_back(PPDI(PI32(l, r), DI(it->second, 1)));
	}

	junctions.assign(mj.begin(), mj.end());
	sbounds.assign(ms.begin(), ms.end());
	tbounds.assign(mt.begin(), mt.end());
	phase.assign(mp.begin(), mp.end());
	reads.assign(mr.begin(), mr.end());

	return 0;
}

int combined_graph::combine_regions(split_interval_map &imap, const combined_graph &gt)
{
	for(int i = 0; i < gt.regions.size(); i++)
	{
		PI32 p = gt.regions[i].first;
		int w = (int)(gt.regions[i].second.first);
		imap += make_pair(ROI(p.first, p.second), w);
	}
	return 0;
}

int combined_graph::combine_junctions(map<PI32, DI> &m, const combined_graph &gt)
{
	for(int i = 0; i < gt.junctions.size(); i++)
	{
		PI32 p = gt.junctions[i].first;
		DI d = gt.junctions[i].second;

		map<PI32, DI>::iterator x = m.find(p);

		if(x == m.end())
		{
			m.insert(pair<PI32, DI>(p, d));
		}
		else 
		{
			x->second.first += d.first;
			x->second.second += d.second;
		}
	}
	return 0;
}

int combined_graph::combine_phase(MV32 &m, const combined_graph &gt)
{
	for(int i = 0; i < gt.phase.size(); i++)
	{
		const vector<int32_t> &p = gt.phase[i].first;
		const vector<int32_t> &v = gt.phase[i].second;
		MV32::iterator x = m.find(p);

		if(x == m.end())
		{
			m.insert(PV32(p, v));
		}
		else 
		{
			x->second.insert(x->second.end(), v.begin(), v.end());
		}
	}
	return 0;
}

int combined_graph::combine_reads(MV32 &m, const combined_graph &gt)
{
	for(int i = 0; i < gt.reads.size(); i++)
	{
		const vector<int32_t> &p = gt.reads[i].first;
		const vector<int32_t> &v = gt.reads[i].second;
		MV32::iterator x = m.find(p);

		if(x == m.end())
		{
			m.insert(PV32(p, v));
		}
		else 
		{
			x->second.insert(x->second.end(), v.begin(), v.end());
		}
	}
	return 0;
}

int combined_graph::combine_start_bounds(map<int32_t, DI> &m, const combined_graph &gt)
{
	for(int i = 0; i < gt.sbounds.size(); i++)
	{
		int32_t p = gt.sbounds[i].first;
		DI d = gt.sbounds[i].second;

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

int combined_graph::combine_end_bounds(map<int32_t, DI> &m, const combined_graph &gt)
{
	for(int i = 0; i < gt.tbounds.size(); i++)
	{
		int32_t p = gt.tbounds[i].first;
		DI d = gt.tbounds[i].second;

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

int combined_graph::resolve(splice_graph &gr, hyper_set &hs, vector<fcluster> &ub)
{
	build_region_index();
	group_junctions();
	build_splice_graph(gr);
	group_start_boundaries(gr);
	group_end_boundaries(gr);
	group_phasing_paths();
	build_phasing_paths(gr, hs);
	refine_splice_graph(gr);
	return 0;
}

int combined_graph::build_region_index()
{
	lindex.clear();
	rindex.clear();
	for(int i = 0; i < regions.size(); i++)
	{
		PI32 p = regions[i].first;
		lindex.insert(pair<int32_t, int>(p.first, i));
		rindex.insert(pair<int32_t, int>(p.second, i));
	}
	return 0;
}

int combined_graph::group_start_boundaries(splice_graph &xr)
{
	smap.clear();
	vector<int> v;
	edge_iterator it;
	PEEI pei = xr.out_edges(0);
	for(edge_iterator it = pei.first; it != pei.second; it++)
	{
		edge_descriptor e = *it;
		assert(e->source() == 0);
		v.push_back(e->target());
	}

	if(v.size() <= 1) return 0;

	sort(v.begin(), v.end());

	int32_t p1 = xr.get_vertex_info(v[0]).lpos;
	int32_t p2 = p1;
	int k1 = v[0];
	int k2 = k1;
	PEB pa = xr.edge(0, v[0]);
	assert(pa.second == true);
	double wa = xr.get_edge_weight(pa.first);
	edge_info ea = xr.get_edge_info(pa.first);

	for(int i = 1; i < v.size(); i++)
	{
		int32_t p = xr.get_vertex_info(v[i]).lpos;
		PEB pb = xr.edge(0, v[i]);
		assert(pb.second == true);
		double wb = xr.get_edge_weight(pb.first);
		edge_info eb = xr.get_edge_info(pb.first);

		bool b = check_continue_vertices(xr, k2, v[i]);

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
				PEB pc = xr.edge(j, j + 1);
				assert(pc.second == true);
				double vc = xr.get_vertex_weight(j);
				double wc = xr.get_edge_weight(pc.first);
				xr.set_vertex_weight(j, vc + wb);
				edge_info ec = xr.get_edge_info(pc.first);
				ec.count += eb.count;
				ec.weight += eb.weight;
				xr.set_edge_weight(pc.first, wc + wb);
				xr.set_edge_info(pc.first, ec);
			}
			wa += wb;
			ea.count += eb.count;
			ea.weight += eb.weight;
			xr.set_edge_weight(pa.first, wa);
			xr.set_edge_info(pa.first, ea);
			xr.remove_edge(pb.first);

			k2 = v[i];
			p2 = p;

			if(meta_verbose >= 2) printf("map start boundary %d:%d (weight = %.2lf) to %d:%d (weight = %.2lf)\n", v[i], p, wb, k1, p1, wa);
		}
	}
	return 0;
}

int combined_graph::group_end_boundaries(splice_graph &xr)
{
	tmap.clear();
	vector<int> v;
	edge_iterator it;
	PEEI pei = xr.in_edges(xr.num_vertices() - 1);
	for(edge_iterator it = pei.first; it != pei.second; it++)
	{
		edge_descriptor e = *it;
		assert(e->target() == xr.num_vertices() - 1);
		v.push_back(e->source());
	}

	if(v.size() <= 1) return 0;

	sort(v.begin(), v.end(), greater<int>());

	int32_t p1 = xr.get_vertex_info(v[0]).rpos;
	int32_t p2 = p1;
	int k1 = v[0];
	int k2 = k1;
	PEB pa = xr.edge(v[0], xr.num_vertices() - 1);
	assert(pa.second == true);
	double wa = xr.get_edge_weight(pa.first);

	for(int i = 1; i < v.size(); i++)
	{
		int32_t p = xr.get_vertex_info(v[i]).rpos;
		PEB pb = xr.edge(v[i], xr.num_vertices() - 1);
		assert(pb.second == true);
		double wb = xr.get_edge_weight(pb.first);

		bool b = check_continue_vertices(xr, v[i], k2);

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
				PEB pc = xr.edge(j, j + 1);
				assert(pc.second == true);
				double wc = xr.get_edge_weight(pc.first);
				double vc = xr.get_vertex_weight(j + 1);
				xr.set_edge_weight(pc.first, wc + wb);
				xr.set_vertex_weight(j + 1, wc + wb);
			}
			wa += wb;
			xr.set_edge_weight(pa.first, wa);
			xr.remove_edge(pb.first);

			k2 = v[i];
			p2 = p;

			if(meta_verbose >= 2) printf("map end boundary %d:%d (weight = %.2lf) to %d:%d (weight = %.2lf)\n", v[i], p, wb, k1, p1, wa);
		}
	}
	return 0;
}

int combined_graph::group_phasing_paths()
{
	for(int i = 0; i < phase.size(); i++)
	{
		vector<int32_t> &v = phase[i].first;
		vector<int32_t> &z = phase[i].second;
		assert(v.size() % 2 == 0);
		assert(z.size() % 3 == 0);

		for(int j = 0; j < z.size(); j++)
		{
			int32_t p = z[j];
			if(j % 3 == 0 && smap.find(p) != smap.end()) z[j] = smap[p];
			if(j % 3 == 1 && tmap.find(p) != tmap.end()) z[j] = tmap[p];
		}

		// TODO merge identical ones to improve speed / save space
	}
	return 0;
}

int combined_graph::group_junctions()
{
	set<int> fb;
	for(int i = 0; i < junctions.size(); i++)
	{
		if(fb.find(i) != fb.end()) continue;
		PPDI x = junctions[i];
		for(int j = i + 1; j < junctions.size(); j++)
		{
			if(fb.find(j) != fb.end()) continue;
			PPDI y = junctions[j];
			double d1 = fabs(x.first.first - y.first.first);
			double d2 = fabs(x.first.second - y.first.second);
			if(d1 + d2 >= max_group_junction_distance) continue;
			if(10 * x.second.first < y.second.first && x.second.second < y.second.second && x.second.second <= 2 && y.second.first <= 100)
			{
				fb.insert(i);
				if(meta_verbose >= 2) printf("filter junction: (%d, %d), w = %.1lf, c = %d -> (%d, %d), w = %.1lf, c = %d\n", 
						x.first.first, x.first.second, x.second.first, x.second.second,
						y.first.first, y.first.second, y.second.first, y.second.second);
			}
			if(x.second.first > 10 * y.second.first && x.second.second > y.second.second && y.second.second <= 2 && y.second.first <= 100)
			{
				if(meta_verbose >= 2) printf("filter junction: (%d, %d), w = %.1lf, c = %d -> (%d, %d), w = %.1lf, c = %d\n",
						y.first.first, y.first.second, y.second.first, y.second.second,
						x.first.first, x.first.second, x.second.first, x.second.second);
				fb.insert(j);
			}
		}
	}

	vector<PPDI> v;
	for(int k = 0; k < junctions.size(); k++)
	{
		PPDI p = junctions[k];
		if(fb.find(k) == fb.end()) v.push_back(p);
	}

	junctions = v;
	return 0;
}

set<int32_t> combined_graph::get_reliable_adjacencies(int samples, double weight)
{
	set<int32_t> s;
	if(regions.size() <= 1) return s;
	for(int i = 0; i < regions.size() - 1; i++)
	{
		int32_t p1 = regions[i + 0].first.second;
		int32_t p2 = regions[i + 1].first.first;
		if(p1 != p2) continue;

		double w1 = regions[i + 0].second.first;
		double w2 = regions[i + 1].second.first;
		int c1 = regions[i + 0].second.second;
		int c2 = regions[i + 1].second.second;

		bool b = false;
		if(w1 >= weight && w2 >= weight) b = true;
		if(c1 >= samples && c2 >= samples) b = true;

		if(b == false) continue;
		s.insert(p1);
	}
	return s;
}

set<int32_t> combined_graph::get_reliable_splices(int samples, double weight)
{
	map<int32_t, DI> m;
	for(int i = 0; i < junctions.size(); i++)
	{
		PI32 p = junctions[i].first;
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

set<PI32> combined_graph::get_reliable_junctions(int samples, double weight)
{
	set<PI32> s;
	for(int i = 0; i < junctions.size(); i++)
	{
		PI32 p = junctions[i].first;
		double w = junctions[i].second.first;
		int c = junctions[i].second.second;

		if(c < samples && w < weight) continue;
		s.insert(p);
	}
	return s;
}

set<int32_t> combined_graph::get_reliable_start_boundaries(int samples, double weight)
{
	map<int32_t, DI> m;
	for(int i = 0; i < sbounds.size(); i++)
	{
		int32_t p = sbounds[i].first;
		int32_t q = p;
		if(smap.find(p) != smap.end()) q = smap[p];
		double w = sbounds[i].second.first;
		int c = sbounds[i].second.second;

		if(m.find(q) == m.end())
		{
			DI di(w, c);
			m.insert(pair<int32_t, DI>(q, di));
		}
		else
		{
			m[q].first += w;
			m[q].second += c;
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

	set<int32_t> ss;
	for(int i = 0; i < sbounds.size(); i++)
	{
		int32_t p = sbounds[i].first;
		int32_t q = p;
		if(smap.find(p) != smap.end()) q = smap[p];
		if(s.find(q) == s.end()) continue;
		ss.insert(p);
	}

	return ss;
}

set<int32_t> combined_graph::get_reliable_end_boundaries(int samples, double weight)
{
	map<int32_t, DI> m;
	for(int i = 0; i < tbounds.size(); i++)
	{
		int32_t p = tbounds[i].first;
		int32_t q = p;
		if(tmap.find(p) != tmap.end()) q = tmap[p];
		double w = tbounds[i].second.first;
		int c = tbounds[i].second.second;

		if(m.find(q) == m.end())
		{
			DI di(w, c);
			m.insert(pair<int32_t, DI>(q, di));
		}
		else
		{
			m[q].first += w;
			m[q].second += c;
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

	set<int32_t> ss;
	for(int i = 0; i < tbounds.size(); i++)
	{
		int32_t p = tbounds[i].first;
		int32_t q = p;
		if(tmap.find(p) != tmap.end()) q = tmap[p];
		if(s.find(q) == s.end()) continue;
		ss.insert(p);
	}

	return ss;
}

int combined_graph::build_splice_graph(splice_graph &xr)
{
	xr.clear();

	xr.gid = gid;
	xr.chrm = chrm;
	xr.strand = strand;

	// add vertices
	xr.add_vertex();	// 0
	PIDI sb = get_leftmost_bound();
	vertex_info v0;
	v0.lpos = sb.first;
	v0.rpos = sb.first;
	xr.set_vertex_weight(0, 0);
	xr.set_vertex_info(0, v0);

	for(int i = 0; i < regions.size(); i++) 
	{
		xr.add_vertex();
		vertex_info vi;
		vi.lpos = regions[i].first.first;
		vi.rpos = regions[i].first.second;
		vi.count = regions[i].second.second;
		double w = regions[i].second.first;
		vi.length = vi.rpos - vi.lpos;
		xr.set_vertex_weight(i + 1, w);
		xr.set_vertex_info(i + 1, vi);
	}

	xr.add_vertex();	// n
	PIDI tb = get_rightmost_bound();
	vertex_info vn;
	vn.lpos = tb.first;
	vn.rpos = tb.first;
	xr.set_vertex_info(regions.size() + 1, vn);
	xr.set_vertex_weight(regions.size() + 1, 0);

	// add sbounds
	for(int i = 0; i < sbounds.size(); i++)
	{
		int32_t p = sbounds[i].first;
		double w = sbounds[i].second.first;
		int c = sbounds[i].second.second;

		assert(lindex.find(p) != lindex.end());
		int k = lindex[p];
		edge_descriptor e = xr.add_edge(0, k + 1);
		edge_info ei;
		ei.weight = w;
		ei.count = c;
		xr.set_edge_info(e, ei);
		xr.set_edge_weight(e, w);
	}

	// add tbounds
	for(int i = 0; i < tbounds.size(); i++)
	{
		int32_t p = tbounds[i].first;
		double w = tbounds[i].second.first;
		int c = tbounds[i].second.second;

		assert(rindex.find(p) != rindex.end());
		int k = rindex[p];
		edge_descriptor e = xr.add_edge(k + 1, regions.size() + 1);
		edge_info ei;
		ei.weight = w;
		ei.count = c;
		xr.set_edge_info(e, ei);
		xr.set_edge_weight(e, w);
	}

	// add junctions
	for(int i = 0; i < junctions.size(); i++)
	{
		PI32 p = junctions[i].first;
		double w = junctions[i].second.first;
		int c = junctions[i].second.second;

		// filtering later on
		//if(parent && c < min_supporting_samples && w < min_splicing_count - 0.01) continue;

		assert(rindex.find(p.first) != rindex.end());
		assert(lindex.find(p.second) != lindex.end());
		int s = rindex[p.first];
		int t = lindex[p.second];
		edge_descriptor e = xr.add_edge(s + 1, t + 1);
		edge_info ei;
		ei.weight = w;
		ei.count = c;
		xr.set_edge_info(e, ei);
		xr.set_edge_weight(e, w);
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
		double w1 = xr.get_out_weights(i + 0);
		double w2 = xr.get_in_weights(i + 1);
		double ws = ss.second.first - w1;
		double wt = tt.second.first - w2;
		double w = (ws + wt) * 0.5;
		*/

		int xd = xr.out_degree(i + 0);
		int yd = xr.in_degree(i + 1);
		double w = (xd < yd) ? ss.second.first : tt.second.first;
		int c = ss.second.second;
		if(ss.second.second > tt.second.second) c = tt.second.second;

		if(w < 1) w = 1;
		edge_descriptor e = xr.add_edge(i + 0, i + 1);
		edge_info ei;
		ei.weight = w;
		ei.count = c;
		xr.set_edge_info(e, ei);
		xr.set_edge_weight(e, w);
	}
	return 0;
}

int combined_graph::build_phasing_paths(splice_graph &gr, hyper_set &hs)
{
	hs.clear();
	for(int i = 0; i < phase.size(); i++)
	{
		vector<int> uu = fetch_vertices_from_coordinates(gr, phase[i].first);
		const vector<int32_t> &z = phase[i].second;
		assert(z.size() % 3 == 0);

		printf("phase %d: core = ( ", i);
		printv(phase[i].first);
		printf(") => ( ");
		printv(uu);
		printf(")\n");

		for(int j = 0; j < z.size() / 3; j++)
		{
			int32_t p1 = z[j * 3 + 0];
			int32_t p2 = z[j * 3 + 1];
			int w = z[j * 3 + 2];

			assert(p1 >= 0 && p2 >= 0);
			assert(lindex.find(p1) != lindex.end());
			assert(rindex.find(p2) != rindex.end());
			int a = lindex[p1] + 1;
			int b = rindex[p2] + 1;

			vector<int> vv;
			if(uu.size() == 0)
			{
				for(int k = a; k <= b; k++) vv.push_back(k);
			}
			else
			{
				for(int k = a; k < uu.front(); k++) vv.push_back(k);
				vv.insert(vv.end(), uu.begin(), uu.end());
				for(int k = uu.back() + 1; k <= b; k++) vv.push_back(k);
			}

			for(int k = 0; k < vv.size(); k++) vv[k]--;

			printf(" (%d, %d, %d) => a = %d, b = %d, ( ", p1, p2, w, a, b);
			printv(vv);
			printf(")\n");
			hs.add_node_list(vv, w);
		}
	}
	return 0;
}

vector<int> combined_graph::fetch_vertices_from_coordinates(splice_graph &gr, const vector<int32_t> &v)
{
	// assume v encodes intron chain coordinates
	vector<int> vv;
	assert(v.size() % 2 == 0);
	if(v.size() <= 0) return vv;

	int n = v.size() / 2;
	vector<PI> pp(n);
	for(int k = 0; k < n; k++)
	{
		int32_t p = v[2 * k + 0];
		int32_t q = v[2 * k + 1];
		assert(p >= 0 && q >= 0);
		assert(p <= q);

		assert(rindex.find(p) != rindex.end());
		assert(lindex.find(q) != lindex.end());
		int kp = rindex[p] + 1;
		int kq = lindex[q] + 1;
		pp[k].first = kp;
		pp[k].second = kq;
	}

	vv.push_back(pp.front().first);
	for(int k = 0; k < n - 1; k++)
	{
		int a = pp[k + 0].second;
		int b = pp[k + 1].first;
		assert(a <= b);
		assert(check_continue_vertices(gr, a, b));
		for(int j = a; j <= b; j++) vv.push_back(j);
	}
	vv.push_back(pp.back().second);

	return vv;
}

int combined_graph::clear()
{
	num_combined = 0;
	gid = "";
	chrm = "";
	strand = '.';
	splices.clear();
	regions.clear();
	junctions.clear();
	sbounds.clear();
	tbounds.clear();
	phase.clear();
	reads.clear();
	lindex.clear();
	rindex.clear();
	smap.clear();
	tmap.clear();
	return 0;
}

int combined_graph::print(int index)
{
	printf("combined-graph %d: #combined = %d, chrm = %s, strand = %c, #regions = %lu, #sbounds = %lu, #tbounds = %lu, #junctions = %lu, #phasing-phase = %lu\n", 
			index, num_combined, chrm.c_str(), strand, regions.size(), sbounds.size(), tbounds.size(), junctions.size(), phase.size());

	//printf("combined-graph %d: #combined = %d, chrm = %s, strand = %c, #regions = %lu, #sbounds = %lu, #tbounds = %lu, #junctions = %lu, #phase = %lu, #reads = %lu\n", 
	//		index, num_combined, chrm.c_str(), strand, regions.size(), sbounds.size(), tbounds.size(), junctions.size(), phase.size(), reads.size());


	for(int i = 0; i < regions.size(); i++)
	{
		PI32 p = regions[i].first;
		DI d = regions[i].second;
		printf("region %d: [%d, %d), w = %.2lf, c = %d\n", i, p.first, p.second, d.first, d.second);
	}
	for(int i = 0; i < junctions.size(); i++)
	{
		PI32 p = junctions[i].first;
		DI d = junctions[i].second;
		printf("junction %d: [%d, %d), w = %.2lf, c = %d, regions = %d -> %d\n", i, p.first, p.second, d.first, d.second, rindex[p.first], lindex[p.second]);
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
	for(int i = 0; i < phase.size(); i++)
	{
		vector<int32_t> &v = phase[i].first;
		vector<int32_t> &z = phase[i].second;
		printf("path %d: %lu components, core = ( ", i, z.size() / 3);
		printv(v);
		printf(")\n");
		for(int k = 0; k < z.size() / 3; k++)
		{
			printf(" bounds = (%d, %d), w = %.2lf, c = 1\n", z[k * 3 + 0], z[k * 3 + 1], 1.0 * z[k * 3 + 2]);
		}
	}
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
