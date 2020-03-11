#include "merged_graph.h"
#include "graph_revise.h"
#include "meta_config.h"
#include "config.h"
#include <sstream>
#include <algorithm>

merged_graph::merged_graph()
{
	num_combined = 0;
}

int merged_graph::solve()
{
	build_region_index();
	group_junctions();

	build_splice_graph(gr);
	group_start_boundaries(gr);
	group_end_boundaries(gr);
	group_phasing_paths();

	refine_splice_graph(gr);
	build_phasing_paths();
	return 0;
}

int merged_graph::clear()
{
	num_combined = 0;
	gid = "";
	chrm = "";
	strand = '.';
	regions.clear();
	junctions.clear();
	sbounds.clear();
	tbounds.clear();
	paths.clear();
	lindex.clear();
	rindex.clear();
	smap.clear();
	tmap.clear();
	gr.clear();
	hs.clear();
	hx.clear();
	return 0;
}

int merged_graph::build_region_index()
{
	for(int i = 0; i < regions.size(); i++)
	{
		PI32 p = regions[i].first;
		lindex.insert(pair<int32_t, int>(p.first, i));
		rindex.insert(pair<int32_t, int>(p.second, i));
	}
	return 0;
}

int merged_graph::group_start_boundaries(splice_graph &xr)
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

		bool b = continue_vertices(k2, v[i], xr);

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

int merged_graph::group_end_boundaries(splice_graph &xr)
{
	tmap.clear();
	vector<int> v;
	edge_iterator it;
	PEEI pei = xr.in_edges(gr.num_vertices() - 1);
	for(edge_iterator it = pei.first; it != pei.second; it++)
	{
		edge_descriptor e = *it;
		assert(e->target() == gr.num_vertices() - 1);
		v.push_back(e->source());
	}

	if(v.size() <= 1) return 0;

	sort(v.begin(), v.end(), greater<int>());

	int32_t p1 = xr.get_vertex_info(v[0]).rpos;
	int32_t p2 = p1;
	int k1 = v[0];
	int k2 = k1;
	PEB pa = xr.edge(v[0], gr.num_vertices() - 1);
	assert(pa.second == true);
	double wa = xr.get_edge_weight(pa.first);

	for(int i = 1; i < v.size(); i++)
	{
		int32_t p = xr.get_vertex_info(v[i]).rpos;
		PEB pb = xr.edge(v[i], gr.num_vertices() - 1);
		assert(pb.second == true);
		double wb = xr.get_edge_weight(pb.first);

		bool b = continue_vertices(v[i], k2, xr);

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

int merged_graph::group_junctions()
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

set<int32_t> merged_graph::get_reliable_adjacencies(int samples, double weight)
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

set<int32_t> merged_graph::get_reliable_splices(int samples, double weight)
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

set<PI32> merged_graph::get_reliable_junctions(int samples, double weight)
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

set<int32_t> merged_graph::get_reliable_start_boundaries(int samples, double weight)
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

set<int32_t> merged_graph::get_reliable_end_boundaries(int samples, double weight)
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

int merged_graph::build_splice_graph(splice_graph &xr)
{
	xr.clear();

	gr.gid = gid;
	gr.chrm = chrm;
	gr.strand = strand;

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

int merged_graph::group_phasing_paths()
{
	for(int i = 0; i < paths.size(); i++)
	{
		vector<int32_t> &v = paths[i].first;
		assert(v.size() >= 2);
		assert(v.size() % 2 == 0);

		int n = v.size() / 2 - 1;

		if(v[0] == -1 && v.size() >= 4)
		{
			assert(v[1] == -1);
			int32_t p = v[2];
			if(smap.find(p) != smap.end()) v[2] = smap[p];
		}
		else if(v[0] > 0)
		{
			int32_t p = v[0];
			if(smap.find(p) != smap.end()) v[0] = smap[p];
		}

		if(v[n * 2] == -2 && v.size() >= 4)
		{
			assert(v[2 * n + 1] == -2);
			int32_t p = v[2 * n - 1];
			if(tmap.find(p) != tmap.end()) v[2 * n - 1] = tmap[p];
		}
		else if(v[n * 2 + 1] > 0)
		{
			int32_t p = v[2 * n + 1];
			if(tmap.find(p) != tmap.end()) v[2 * n + 1] = tmap[p];
		}

	}

	if(paths.size() <= 1) return 0;

	sort(paths.begin(), paths.end(), compare_phasing_path);

	vector<PVDI> vv;
	PVDI p = paths[0];

	for(int i = 1; i < paths.size(); i++)
	{
		if(paths[i].first == p.first)
		{
			p.second.first += paths[i].second.first;
			p.second.second += paths[i].second.second;
		}
		else
		{
			vv.push_back(p);
			p = paths[i];
		}
	}

	paths = vv;
	return 0;
}

int merged_graph::build_phasing_paths()
{
	hs.clear();
	for(int i = 0; i < paths.size(); i++)
	{
		const vector<int32_t> &v = paths[i].first;
		assert(v.size() >= 2);
		assert(v.size() % 2 == 0);

		double w = paths[i].second.first;
		int c = paths[i].second.second;

		// filtering here
		//if((v.front() != -1 || v.back() != -2) &&)
		//if(c < min_supporting_samples &&);
		//if(w < min_phasing_count) continue;

		vector<int> vv;
		int n = v.size() / 2 - 1;
		for(int k = 0; k <= n; k++)
		{
			if(v[2 * k + 0] == -1)
			{
				assert(k == 0);
				assert(v[2 * k + 1] == -1);
				vv.push_back(0); 
			}
			else if(v[2 * k + 0] == -2)
			{
				assert(v[2 * k + 1] == -2);
				vv.push_back(gr.num_vertices() - 1); 
			}
			else
			{
				assert(v[2 * k + 0] >= 0);
				assert(v[2 * k + 1] >= 0);
				int32_t p = v[2 * k + 0];
				int32_t q = v[2 * k + 1];

				assert(lindex.find(p) != lindex.end());
				assert(rindex.find(q) != rindex.end());
				int kp = lindex[p] + 1;
				int kq = rindex[q] + 1;
				if(vv.size() >= 1) assert(kp > vv.back());
				for(int j = kp; j <= kq; j++)
				{
					vv.push_back(j);
					if(j == kp) continue;
					assert(regions[j - 2].first.second == regions[j - 1].first.first);
				}
			}
		}
		
		for(int k = 0; k < vv.size(); k++) vv[k]--;

		if(meta_verbose >= 2)
		{
			printf("phasing path %d: count = %d, weight = %.2lf, vertices = ", i, c, w);
			printv(vv);
			printf("\n");
		}
		hs.add_node_list(vv, (int)w);
	}
	return 0;
}


int merged_graph::build_phasing_paths(const vector<PVDI> &px)
{
	hx.clear();

	int32_t leftbound = get_leftmost_bound().first;
	int32_t rightbound = get_rightmost_bound().first;

	for(int i = 0; i < px.size(); i++)
	{
		const vector<int32_t> &v = px[i].first;
		assert(v.size() >= 2);
		assert(v.size() % 2 == 0);

		/*
		if(v.front() < 0 && v.back() < 0 && v.size() <= 6) continue;
		if(v.front() < 0 && v.back() >= 0 && v.size() <= 4) continue;
		if(v.front() >=0 && v.back() < 0 && v.size() <= 4) continue;
		if(v.front() >=0 && v.back() >=0 && v.size() <= 2) continue;
		*/

		if(v.front() >= 0 && v.front() > rightbound) continue;
		if(v.back() >= 0 && v.back() < leftbound) continue;
		if(v.front() < 0 && v[2] > rightbound) continue;
		if(v.back() < 0 && v[v.size() - 3] < leftbound) continue;

		double w = px[i].second.first;
		int c = px[i].second.second;

		vector<int> vv;
		int n = v.size() / 2 - 1;
		bool fail = false;
		for(int k = 0; k <= n; k++)
		{
			if(v[2 * k + 0] == -1)
			{
				assert(k == 0);
				assert(v[2 * k + 1] == -1);
				vv.push_back(0);
			}
			else if(v[2 * k + 0] == -2)
			{
				assert(v[2 * k + 1] == -2);
				vv.push_back(gr.num_vertices() - 1);
			}
			else if((k == 0 && v[0] != -1) || (k == 1 && v[0] == -1))
			{
				assert(v[2 * k + 0] >= 0);
				assert(v[2 * k + 1] >= 0);
				int32_t p = v[2 * k + 0];
				int32_t q = v[2 * k + 1];

				if(rindex.find(q) == rindex.end()) fail = true;
				if(fail == true) break;
				int kq = rindex[q] + 1;
				//int kp = kq;

				int r = locate_left_region(p, 0, regions.size());
				if(r == -1) fail = true;
				if(fail == true) break;
				int kp = r + 1;

				if(vv.size() >= 1) assert(kp > vv.back());
				for(int j = kp; j <= kq; j++)
				{
					vv.push_back(j);
				}
			}
			else if((k == n && v[2 * n] != -2) || (k == n - 1 && v[2 * n] == -2))
			{
				assert(v[2 * k + 0] >= 0);
				assert(v[2 * k + 1] >= 0);
				int32_t p = v[2 * k + 0];
				int32_t q = v[2 * k + 1];

				if(lindex.find(p) == lindex.end()) fail = true;
				if(fail == true) break;
				int kp = lindex[p] + 1;
				//int kq = kp;

				int r = locate_right_region(q, 0, regions.size());
				if(r == -1) fail = true;
				if(fail == true) break;
				int kq = r + 1;

				if(vv.size() >= 1) assert(kp > vv.back());
				for(int j = kp; j <= kq; j++)
				{
					vv.push_back(j);
				}
			}
			else
			{
				assert(v[2 * k + 0] >= 0);
				assert(v[2 * k + 1] >= 0);
				int32_t p = v[2 * k + 0];
				int32_t q = v[2 * k + 1];

				if(lindex.find(p) == lindex.end()) fail = true;
				if(fail == true) break;

				if(rindex.find(q) == rindex.end()) fail = true;
				if(fail == true) break;

				int kp = lindex[p] + 1;
				int kq = rindex[q] + 1;
				if(vv.size() >= 1) assert(kp > vv.back());
				for(int j = kp; j <= kq; j++)
				{
					vv.push_back(j);
				}
			}
		}

		if(fail == true) continue;
		
		for(int k = 0; k < vv.size(); k++) vv[k]--;

		if(meta_verbose >= 2)
		{
			printf("phasing path %d: count = %d, weight = %.2lf, vertices = ", i, c, w);
			printv(vv);
			printf("\n");
		}

		hx.add_node_list(vv, (int)w);
	}
	return 0;
}

int merged_graph::build(combined_graph &cb)
{
	gid = "TODO";
	chrm = cb.chrm;
	strand = cb.strand;
	num_combined = cb.num_combined;

	for(SIMI it = cb.imap.begin(); it != cb.imap.end(); it++)
	{
		int32_t l = lower(it->first);
		int32_t r = upper(it->first);
		regions.push_back(PPDI(PI32(l, r), DI(it->second, 1)));
	}

	for(map<int32_t, DI>::iterator it = cb.sbounds.begin(); it != cb.sbounds.end(); it++)
	{
		int32_t p = it->first;
		double w = it->second.first;
		int c = it->second.second;
		sbounds.push_back(PIDI(p, DI(w, c)));
	}

	for(map<int32_t, DI>::iterator it = cb.tbounds.begin(); it != cb.tbounds.end(); it++)
	{
		int32_t p = it->first;
		double w = it->second.first;
		int c = it->second.second;
		tbounds.push_back(PIDI(p, DI(w, c)));
	}

	for(map<PI32, DI>::iterator it = cb.junctions.begin(); it != cb.junctions.end(); it++)
	{
		int32_t s = it->first.first;
		int32_t t = it->first.second;
		double w = it->second.first;
		int c = it->second.second;

		if(s >= t) continue;
		junctions.push_back(PPDI(PI32(s, t), DI(w, c)));
	}

	for(map<vector<int32_t>, DI>::iterator it = cb.phase.begin(); it != cb.phase.end(); it++)
	{
		const vector<int32_t> &vv = it->first;
		double w = it->second.first;
		int c = it->second.second;
		paths.push_back(PVDI(vv, DI(w, c)));
	}

	for(map<vector<int32_t>, DI>::iterator it = cb.paths.begin(); it != cb.paths.end(); it++)
	{
		const vector<int32_t> &vv = it->first;
		double w = it->second.first;
		int c = it->second.second;
		paths.push_back(PVDI(vv, DI(w, c)));
	}

	return 0;
}

int merged_graph::build(istream &is, const string &id, const string &ch, char st, int num)
{
	gid = id;
	chrm = ch;
	strand = st;
	num_combined = num;

	char line[10240];
	char name[10240];
	while(is.getline(line, 10240, '\n'))
	{
		stringstream sstr(line);
		if(string(line).length() == 0) break;
		
		sstr >> name;
		if(string(name) == "region")
		{
			double weight;
			int32_t lpos;
			int32_t rpos;
			sstr >> lpos >> rpos >> weight;
			regions.push_back(PPDI(PI32(lpos, rpos), DI(weight, 1)));
		}
		else if(string(name) == "sbound")
		{
			int32_t p;
			double w;
			int c;
			sstr >> p >> w >> c;
			sbounds.push_back(PIDI(p, DI(w, c)));
		}
		else if(string(name) == "tbound")
		{
			int32_t p;
			double w;
			int c;
			sstr >> p >> w >> c;
			tbounds.push_back(PIDI(p, DI(w, c)));
		}
		else if(string(name) == "junction")
		{
			int32_t s, t;
			double w;
			int c;
			sstr >> s >> t >> w >> c;
			junctions.push_back(PPDI(PI32(s, t), DI(w, c)));
		}
		else if(string(name) == "phase")
		{
			int z;
			sstr >> z;
			vector<int32_t> v;
			v.resize(z);
			for(int k = 0; k < z; k++) sstr >> v[k];
			double w;
			int c;
			sstr >> w >> c;
			paths.push_back(PVDI(v, DI(w, c)));
		}
		else if(string(name) == "path" && parent == true)
		{
			int z;
			sstr >> z;
			vector<int32_t> v;
			v.resize(z);
			for(int k = 0; k < z; k++) sstr >> v[k];
			double w;
			int c;
			sstr >> w >> c;
			paths.push_back(PVDI(v, DI(w, c)));
		}
		else
		{
			break;
		}
	}
	return 0;
}

int merged_graph::print(int index)
{
	printf("combined-graph %d: #combined = %d, chrm = %s, strand = %c, #regions = %lu, #sbounds = %lu, #tbounds = %lu, #junctions = %lu, #phasing-paths = %lu\n", 
			index, num_combined, chrm.c_str(), strand, regions.size(), sbounds.size(), tbounds.size(), junctions.size(), paths.size());

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
	for(int i = 0; i < paths.size(); i++)
	{
		vector<int32_t> &v = paths[i].first;
		DI d = paths[i].second;
		printf("path %d: w = %.2lf, c = %d, list = ", i, d.first, d.second);
		for(int k = 0; k < v.size(); k++) printf("%d ", v[k]);
		printf("\n");
	}
	return 0;
}


bool merged_graph::continue_vertices(int x, int y, splice_graph &xr)
{
	if(x >= y) return true;
	for(int i = x; i < y; i++)
	{
		PEB p = xr.edge(i, i + 1);
		if(p.second == false) return false;
		if(xr.get_vertex_info(i).rpos != xr.get_vertex_info(i + 1).lpos) return false;
	}
	return true;
}

PIDI merged_graph::get_leftmost_bound()
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

PIDI merged_graph::get_rightmost_bound()
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

int merged_graph::locate_left_region(int32_t p, int kl, int kr)
{
	if(kl >= kr) return -1;

	if(kl + 1 == kr)
	{
		PPDI &r = regions[kl];
		if(p < r.first.second) return kl;
		else return -1;
	}

	int m = (kl + kr) / 2;
	assert(kl < m && m < kr);

	if(p < regions[m - 1].first.second) return locate_left_region(p, kl, m);
	else return locate_left_region(p, m, kr);
}

int merged_graph::locate_right_region(int32_t p, int kl, int kr)
{
	if(kl >= kr) return -1;

	if(kl + 1 == kr)
	{
		PPDI &r = regions[kl];
		if(p >= r.first.first) return kl;
		else return -1;
	}

	int m = (kl + kr) / 2;
	assert(kl < m && m < kr);

	if(p >= regions[m].first.first) return locate_right_region(p, m, kr);
	else return locate_right_region(p, kl, m);
}

bool compare_phasing_path(const PVDI &x, const PVDI &y)
{
	if(x.first < y.first) return true;
	if(x.first > y.first) return false;
	if(x.second.second < y.second.second) return true;
	if(x.second.second > y.second.second) return false;
	return x.second.first < y.second.first;
}
