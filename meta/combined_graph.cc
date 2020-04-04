#include "meta_config.h"
#include "combined_graph.h"
#include "config.h"
#include "essential.h"
#include <sstream>
#include <algorithm>

combined_graph::combined_graph()
{
	num_combined = 0;
	strand = '?';
}

int combined_graph::build(splice_graph &gr, const phase_set &ps, const vector<pereads_cluster> &ub)
{
	chrm = gr.chrm;
	strand = gr.strand;
	num_combined = 1;

	build_regions(gr);
	build_start_bounds(gr);
	build_end_bounds(gr);
	build_splices_junctions(gr);
	phases = ps;
	preads = ub;
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
		if(p1 >= p2) continue;

		PI32 p(p1, p2);
		junctions.push_back(PPDI(p, DI(w, 1)));
		sp.insert(p1);
		sp.insert(p2);
	}
	splices.assign(sp.begin(), sp.end());
	sort(splices.begin(), splices.end());
	return 0;
}

/* TODO TODO
int combined_graph::combine_extra_bridged_reads(const vector< vector<int32_t> > &exon_chains, const vector<int> &weights)
{
	assert(exon_chains.size() == weights.size());
	if(exon_chains.size() == 0) return 0;
	int32_t lb = get_leftmost_bound().first;
	int32_t rb = get_rightmost_bound().first;

	// index phase
	map<vector<int32_t>, int> mm;
	for(int k = 0; k < phase.size(); k++)
	{
		mm.insert(pair<vector<int32_t>, int>(phase[k].vv, k));
	}

	// create a new combined graph to be combined
	combined_graph cb;
	for(int k = 0; k < exon_chains.size(); k++)
	{
		const vector<int32_t> &vv = exon_chains[k];
		int c = weights[k];
		assert(vv.size() % 2 == 0);
		if(vv.size() <= 1) continue;
		if(vv.front() < lb) continue;
		if(vv.back() > rb) continue;

		// regions
		for(int i = 0; i < vv.size() / 2; i++)
		{
			PI32 p(vv[i * 2 + 0], vv[i * 2 + 1]);
			cb.regions.push_back(PPDI(p, DI(c, 1)));
		}

		// junctions
		for(int i = 0; i < vv.size() / 2 - 1; i++)
		{
			PI32 p(vv[i * 2 + 1], vv[i * 2 + 2]);
			cb.junctions.push_back(PPDI(p, DI(c, 1)));
		}

		// phase
		vector<int32_t> zz(vv.begin() + 1, vv.end() - 1);
		map<vector<int32_t>, int>::iterator tp = mm.find(zz);
		if(tp == mm.end()) 
		{
			rcluster r;
			r.vv = zz;
			r.vl.push_back(vv.front());
			r.vr.push_back(vv.back());
			r.cc.push_back(c);
			mm.insert(pair<vector<int32_t>, int>(zz, phase.size()));
			phase.push_back(r);
		}
		else 
		{
			int k = tp->second;
			assert(zz == phase[k].vv);
			phase[k].vl.push_back(vv.front());
			phase[k].vr.push_back(vv.back());
			phase[k].cc.push_back(c);
		}
	}

	split_interval_map imap;
	map<PI32, DI> mj;

	combine_regions(imap, *this);
	combine_regions(imap, cb);
	combine_junctions(mj, *this);
	combine_junctions(mj, cb);

	regions.clear();
	for(SIMI it = imap.begin(); it != imap.end(); it++)
	{
		int32_t l = lower(it->first);
		int32_t r = upper(it->first);
		regions.push_back(PPDI(PI32(l, r), DI(it->second, 1)));
	}
	junctions.assign(mj.begin(), mj.end());

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
*/

int combined_graph::get_overlapped_splice_positions(const vector<int32_t> &v) const
{
	vector<int32_t> vv(v.size(), 0);
	vector<int32_t>::iterator it = set_intersection(v.begin(), v.end(), splices.begin(), splices.end(), vv.begin());
	return it - vv.begin();
}

/*
int combined_graph::combine_children()
{
	if(children.size() == 0) return 0;

	split_interval_map imap;
	map<PI32, DI> mj;
	map<int32_t, DI> ms;
	map<int32_t, DI> mt;
	phases.clear();
	preads.clear();

	int num = 0;
	for(int i = 0; i < children.size(); i++)
	{
		combined_graph &gt = children[i];
		combine_regions(imap, gt);
		combine_junctions(mj, gt);
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

	return 0;
}
*/

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

int combined_graph::combine_junctions(map<PI32, DI> &m) const
{
	for(int i = 0; i < junctions.size(); i++)
	{
		PI32 p = junctions[i].first;
		DI d = junctions[i].second;

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

int combined_graph::build_splice_graph(splice_graph &gr)
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

	// add junctions
	for(int i = 0; i < junctions.size(); i++)
	{
		PI32 p = junctions[i].first;
		double w = junctions[i].second.first;
		int c = junctions[i].second.second;

		assert(gr.rindex.find(p.first) != gr.rindex.end());
		assert(gr.lindex.find(p.second) != gr.lindex.end());
		int s = gr.rindex[p.first];
		int t = gr.lindex[p.second];
		edge_descriptor e = gr.add_edge(s, t);
		edge_info ei;
		ei.weight = w;
		ei.count = c;
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
		int c = ss.second.second;
		if(ss.second.second > tt.second.second) c = tt.second.second;

		if(w < 1) w = 1;
		edge_descriptor e = gr.add_edge(i + 0, i + 1);
		edge_info ei;
		ei.weight = w;
		ei.count = c;
		gr.set_edge_info(e, ei);
		gr.set_edge_weight(e, w);
	}
	return 0;
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
	phases.clear();
	preads.clear();
	return 0;
}

int combined_graph::print(int index)
{
	printf("combined-graph %d: #combined = %d, chrm = %s, strand = %c, #regions = %lu, #sbounds = %lu, #tbounds = %lu, #junctions = %lu, #phasing-phase = %lu\n", 
			index, num_combined, chrm.c_str(), strand, regions.size(), sbounds.size(), tbounds.size(), junctions.size(), phases.pmap.size());

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
		printf("junction %d: [%d, %d), w = %.2lf, c = %d\n", i, p.first, p.second, d.first, d.second);
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
	phases.print();
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
