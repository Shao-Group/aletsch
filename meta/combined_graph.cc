#include "combined_graph.h"
#include "draw.h"

combined_graph::combined_graph()
{
	num_combined = 0;
}

combined_graph::combined_graph(const string &line)
{
	hline = line;
	num_combined = 0;
}

int combined_graph::combine(const combined_graph &gt)
{
	if(children.size() == 0) children.push_back(*this);
	if(gt.children.size() == 0) children.push_back(gt);
	children.insert(children.end(), gt.children.begin(), gt.children.end());

	if(chrm == "") chrm = gt.chrm;
	if(strand == '?') strand = gt.strand;
	assert(gt.chrm == chrm);
	assert(gt.strand == strand);

	combine_regions(gt);
	combine_junctions(gt);
	combine_phase(gt);
	combine_paths(gt);
	combine_start_bounds(gt);
	combine_end_bounds(gt);
	combine_splice_positions(gt);
	num_combined += gt.num_combined;
	return 0;
}

int combined_graph::combine_regions(const combined_graph &gt)
{
	for(SIMI it = gt.imap.begin(); it != gt.imap.end(); it++)
	{
		imap += make_pair(it->first, it->second);
	}
	return 0;
}

int combined_graph::combine_junctions(const combined_graph &gt)
{
	for(map<PI32, DI>::const_iterator it = gt.junctions.begin(); it != gt.junctions.end(); it++)
	{
		PI32 p = it->first;
		DI d = it->second;

		map<PI32, DI>::iterator x = junctions.find(p);

		if(x == junctions.end())
		{
			junctions.insert(pair<PI32, DI>(p, d));
		}
		else 
		{
			x->second.first += d.first;
			x->second.second += d.second;
		}
	}
	return 0;
}

int combined_graph::combine_phase(const combined_graph &gt)
{
	for(map<vector<int32_t>, DI>::const_iterator it = gt.phase.begin(); it != gt.phase.end(); it++)
	{
		const vector<int32_t> &p = it->first;
		DI d = it->second;

		map<vector<int32_t>, DI>::iterator x = phase.find(p);

		if(x == phase.end())
		{
			phase.insert(pair<vector<int32_t>, DI>(p, d));
		}
		else 
		{
			x->second.first += d.first;
			x->second.second += d.second;
		}
	}
	return 0;
}

int combined_graph::combine_paths(const combined_graph &gt)
{
	for(map<vector<int32_t>, DI>::const_iterator it = gt.paths.begin(); it != gt.paths.end(); it++)
	{
		const vector<int32_t> &p = it->first;
		DI d = it->second;

		map<vector<int32_t>, DI>::iterator x = paths.find(p);

		if(x == paths.end())
		{
			paths.insert(pair<vector<int32_t>, DI>(p, d));
		}
		else 
		{
			x->second.first += d.first;
			x->second.second += d.second;
		}
	}
	return 0;
}

int combined_graph::combine_start_bounds(const combined_graph &gt)
{
	for(map<int32_t, DI>::const_iterator it = gt.sbounds.begin(); it != gt.sbounds.end(); it++)
	{
		int32_t p = it->first;
		DI d = it->second;

		map<int32_t, DI>::iterator x = sbounds.find(p);

		if(x == sbounds.end())
		{
			sbounds.insert(pair<int32_t, DI>(p, d));
		}
		else 
		{
			x->second.first += d.first;
			x->second.second += d.second;
		}
	}
	return 0;
}

int combined_graph::combine_end_bounds(const combined_graph &gt)
{
	for(map<int32_t, DI>::const_iterator it = gt.tbounds.begin(); it != gt.tbounds.end(); it++)
	{
		int32_t p = it->first;
		DI d = it->second;

		map<int32_t, DI>::iterator x = tbounds.find(p);

		if(x == tbounds.end())
		{
			tbounds.insert(pair<int32_t, DI>(p, d));
		}
		else 
		{
			x->second.first += d.first;
			x->second.second += d.second;
		}
	}
	return 0;
}

int combined_graph::combine_splice_positions(const combined_graph &gt)
{
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

int combined_graph::build(splice_graph &gr, hyper_set &hs)
{
	chrm = gr.chrm;
	strand = gr.strand;
	int n = gr.num_vertices() - 1;

	// vertices
	for(int i = 1; i < n; i++)
	{
		double weight = gr.get_vertex_weight(i);
		vertex_info vi = gr.get_vertex_info(i);
		imap += make_pair(ROI(vi.lpos, vi.rpos), (int)(weight));
	}

	// sbound
	PEEI pei = gr.out_edges(0);
	for(edge_iterator it = pei.first; it != pei.second; it++)
	{
		int s = (*it)->source(); 
		int t = (*it)->target();
		assert(s == 0 && t > s);
		if(t == n) continue;
		double w = gr.get_edge_weight(*it);
		int32_t p = gr.get_vertex_info(t).lpos;
		int c = 1;

		map<int32_t, DI>::iterator tp = sbounds.find(p);
		if(tp == sbounds.end()) 
		{
			sbounds.insert(pair<int32_t, DI>(p, DI(w, c)));
		}
		else 
		{
			tp->second.first += w;
			tp->second.second += c;
		}
	}

	// tbound
	pei = gr.in_edges(n);
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

		map<int32_t, DI>::iterator tp = tbounds.find(p);
		if(tp == tbounds.end()) 
		{
			tbounds.insert(pair<int32_t, DI>(p, DI(w, c)));
		}
		else 
		{
			tp->second.first += w;
			tp->second.second += c;
		}
	}

	// junctions and splices
	pei = gr.edges();
	set<int32_t> sp;
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
		map<PI32, DI>::iterator tp = junctions.find(p);
		if(tp == junctions.end()) 
		{
			junctions.insert(pair<PI32, DI>(p, DI(w, c)));
		}
		else 
		{
			tp->second.first += w;
			tp->second.second += c;
		}
		sp.insert(p1);
		sp.insert(p2);
	}
	splices.clear();
	splices.assign(sp.begin(), sp.end());
	sort(splices.begin(), splices.end());

	// phasing paths
	for(MVII::const_iterator it = hs.nodes.begin(); it != hs.nodes.end(); it++)
	{
		const vector<int> &v = it->first;
		int w = it->second;
		int c = 1;

		if(v.size() <= 0) continue;
		vector<int32_t> vv;
		int32_t pre = -9999;

		if(v.front() == 0) vv.push_back(-1);
		if(v.front() == 0) vv.push_back(-1);

		for(int i = 0; i < v.size(); i++)
		{
			int p = v[i];

			if(p == 0) continue;
			if(p == n) continue;

			int32_t pp = gr.get_vertex_info(p).lpos;
			int32_t qq = gr.get_vertex_info(p).rpos;

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

		map<vector<int32_t>, DI>::iterator tp = phase.find(vv);

		if(tp == phase.end()) 
		{
			phase.insert(pair<vector<int32_t>, DI>(vv, DI(w, c)));
		}
		else 
		{
			tp->second.first += w;
			tp->second.second += c;
		}
	}

	return 0;
}

int combined_graph::build(istream &is, const string &ch, char st)
{
	chrm = ch;
	strand = st;

	char line[10240];
	char name[10240];

	set<int32_t> sp(splices.begin(), splices.end());
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
			imap += make_pair(ROI(lpos, rpos), (int)(weight));
		}
		else if(string(name) == "sbound")
		{
			int32_t p;
			double w;
			int c;
			sstr >> p >> w >> c;

			map<int32_t, DI>::iterator it = sbounds.find(p);
			if(it == sbounds.end()) 
			{
				sbounds.insert(pair<int32_t, DI>(p, DI(w, c)));
			}
			else 
			{
				it->second.first += w;
				it->second.second += c;
			}
		}
		else if(string(name) == "tbound")
		{
			int32_t p;
			double w;
			int c;
			sstr >> p >> w >> c;
			map<int32_t, DI>::iterator it = tbounds.find(p);
			if(it == tbounds.end()) 
			{
				tbounds.insert(pair<int32_t, DI>(p, DI(w, c)));
			}
			else 
			{
				it->second.first += w;
				it->second.second += c;
			}
		}
		else if(string(name) == "junction")
		{
			int32_t s, t;
			double w;
			int c;
			sstr >> s >> t >> w >> c;

			assert(s < t);
			sp.insert(s);
			sp.insert(t);

			PI32 p(s, t);
			map<PI32, DI>::iterator it = junctions.find(p);

			if(it == junctions.end()) 
			{
				junctions.insert(pair<PI32, DI>(p, DI(w, c)));
			}
			else 
			{
				it->second.first += w;
				it->second.second += c;
			}
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

			map<vector<int32_t>, DI>::iterator it = phase.find(v);

			if(it == phase.end()) 
			{
				phase.insert(pair<vector<int32_t>, DI>(v, DI(w, 1)));
			}
			else 
			{
				it->second.first += w;
				it->second.second += c;
			}
		}
		else if(string(name) == "path" && false)	// TODO
		{
			int z;
			sstr >> z;
			vector<int32_t> v;
			v.resize(z);
			for(int k = 0; k < z; k++) sstr >> v[k];
			double w;
			int c;
			sstr >> w >> c;

			map<vector<int32_t>, DI>::iterator it = paths.find(v);

			if(it == paths.end()) 
			{
				paths.insert(pair<vector<int32_t>, DI>(v, DI(w, 1)));
			}
			else 
			{
				it->second.first += w;
				it->second.second += c;
			}
		}

		else
		{
			break;
		}
	}

	splices.clear();
	splices.assign(sp.begin(), sp.end());
	sort(splices.begin(), splices.end());
	num_combined++;

	return 0;
}

int combined_graph::write(ostream &os)
{
	os<<fixed;
	os.precision(2);

	for(SIMI it = imap.begin(); it != imap.end(); it++)
	{
		int32_t l = lower(it->first);
		int32_t r = upper(it->first);
		os << "region " << l << " " << r << " " << it->second << endl;
	}

	for(map<int32_t, DI>::iterator it = sbounds.begin(); it != sbounds.end(); it++)
	{
		int32_t p = it->first;
		double w = it->second.first;
		int c = it->second.second;
		os << "sbound " << p << " " << w << " " << c << endl;
	}

	for(map<int32_t, DI>::iterator it = tbounds.begin(); it != tbounds.end(); it++)
	{
		int32_t p = it->first;
		double w = it->second.first;
		int c = it->second.second;
		os << "tbound " << p << " " << w << " " << c << endl;
	}

	for(map<PI32, DI>::iterator it = junctions.begin(); it != junctions.end(); it++)
	{
		int32_t s = it->first.first;
		int32_t t = it->first.second;
		double w = it->second.first;
		int c = it->second.second;

		if(s >= t) continue;
		os << "junction " << s << " " << t << " " << w << " " << c << endl;
	}

	for(map<vector<int32_t>, DI>::iterator it = phase.begin(); it != phase.end(); it++)
	{
		const vector<int32_t> &vv = it->first;
		double w = it->second.first;
		int c = it->second.second;

		os << "phase " << vv.size();
		for(int i = 0; i < vv.size(); i++) os << " " << vv[i];
		os << " " << w << " " << c << endl;
	}

	for(map<vector<int32_t>, DI>::iterator it = paths.begin(); it != paths.end(); it++)
	{
		const vector<int32_t> &vv = it->first;
		double w = it->second.first;
		int c = it->second.second;

		os << "path " << vv.size();
		for(int i = 0; i < vv.size(); i++) os << " " << vv[i];
		os << " " << w << " " << c << endl;
	}
	return 0;
}

int combined_graph::write(ostream &os, int index, bool headers)
{
	char name[10240];
	sprintf(name, "graph.%d", index);
	int n = std::distance(imap.begin(), imap.end());
	os << "# " << name << " " << chrm.c_str() << " " << strand << " " << num_combined;
	os << " " << n << " " << sbounds.size() << " " << tbounds.size() << " " << junctions.size() << " " << phase.size() << " " << paths.size() << endl;
	if(headers == true) return 0;

	// write current graph
	write(os);
	os << endl;

	// write children
	for(int k = 0; k < children.size(); k++)
	{
		combined_graph &cb = children[k];
		os << cb.hline.c_str() << endl;
		cb.write(os);
		os << endl;
	}
	return 0;
}

PI32 combined_graph::get_bounds()
{
	if(sbounds.size() == 0 || tbounds.size() == 0) return PI32(-1, -1);
	map<int32_t, DI>::iterator x1 = sbounds.begin();
	map<int32_t, DI>::iterator x2 = tbounds.end();
	x2--;
	return PI32(x1->first, x2->first);
}

int combined_graph::print(int index)
{
	PI32 p = get_bounds();
	printf("combined-graph %d: #combined = %d, children = %lu, chrm = %s, strand = %c, #regions = %lu, #sbounds = %lu, #tbounds = %lu, #junctions = %lu, #phase = %lu, #splices = %lu, $paths = %lu, boundary = [%d, %d)\n", 
			index, num_combined, children.size(), chrm.c_str(), strand, std::distance(imap.begin(), imap.end()), sbounds.size(), tbounds.size(), junctions.size(), phase.size(), splices.size(), paths.size(), p.first, p.second);
	return 0;
}

int combined_graph::analyze(int index)
{
	return 0;
}
