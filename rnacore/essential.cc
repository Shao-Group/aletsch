#include "essential.h"
#include "splice_graph.h"

int build_child_splice_graph(splice_graph &root, splice_graph &gr, const set<int> &ss)
{
	// make sure ss does not contain 0 and n
	gr.clear();
	if(ss.size() <= 0) return 0;

	vector<int> vv(ss.begin(), ss.end());
	sort(vv.begin(), vv.end());
	map<int, int> a2b;
	map<int, int> b2a;
	for(int i = 0; i < vv.size(); i++)
	{
		a2b.insert(pair<int, int>(vv[i], i + 1));
		b2a.insert(pair<int, int>(i + 1, vv[i]));
	}

	int32_t lpos = root.get_vertex_info(vv.front()).lpos;
	int32_t rpos = root.get_vertex_info(vv.back()).rpos;

	// vertices
	gr.add_vertex();
	vertex_info vi0;
	vi0.lpos = lpos;
	vi0.rpos = lpos;
	gr.set_vertex_weight(0, 0);
	gr.set_vertex_info(0, vi0);

	for(int i = 0; i < vv.size(); i++)
	{
		int k = vv[i];
		gr.add_vertex();
		gr.set_vertex_weight(i + 1, root.get_vertex_weight(k));
		gr.set_vertex_info(i + 1, root.get_vertex_info(k));
	}

	gr.add_vertex();
	vertex_info vin;
	vin.lpos = rpos;
	vin.rpos = rpos;
	gr.set_vertex_weight(vv.size() + 1, 0);
	gr.set_vertex_info(vv.size() + 1, vin);

	// edges
	edge_iterator it1, it2;
	PEEI pei;
	for(pei = root.out_edges(0), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		int t = (*it1)->target();
		if(ss.find(t) == ss.end()) continue;
		assert(a2b.find(t) != a2b.end());
		int y = a2b[t];

		edge_descriptor e = gr.add_edge(0, y);
		gr.set_edge_weight(e, root.get_edge_weight(*it1));
		gr.set_edge_info(e, root.get_edge_info(*it1));
	}

	int n = root.num_vertices() - 1;
	for(int i = 0; i < vv.size(); i++)
	{
		int s = vv[i];
		assert(s != 0 && s != n);
		assert(a2b.find(s) != a2b.end());
		int x = a2b[s];

		pei = root.out_edges(s); 
		for(it1 = pei.first; it1 != pei.second; it1++)
		{
			int t = (*it1)->target();
			assert(t == n || ss.find(t) != ss.end());
			assert(t == n || a2b.find(t) != a2b.end());
			int y = ((t == n) ? gr.num_vertices() - 1 : a2b[t]);

			edge_descriptor e = gr.add_edge(x, y);
			gr.set_edge_weight(e, root.get_edge_weight(*it1));
			gr.set_edge_info(e, root.get_edge_info(*it1));
		}
	}
	return 0;
}

int build_child_hyper_set(hyper_set &hyper, hyper_set &hs, const set<int> &ss)
{
	// hyper-set
	hs.clear();
	if(ss.size() <= 0) return 0;

	vector<int> vv(ss.begin(), ss.end());
	sort(vv.begin(), vv.end());
	map<int, int> a2b;
	map<int, int> b2a;
	for(int i = 0; i < vv.size(); i++)
	{
		a2b.insert(pair<int, int>(vv[i], i + 1));
		b2a.insert(pair<int, int>(i + 1, vv[i]));
	}

	for(MVII::const_iterator it = hyper.nodes.begin(); it != hyper.nodes.end(); it++)
	{
		vector<int> v = it->first;
		int c = it->second;

		bool b = true;
		vector<int> vv;
		for(int k = 0; k < v.size(); k++)
		{
			if(ss.find(v[k]) == ss.end()) b = false;
			if(b == false) break;
			assert(a2b.find(v[k]) != a2b.end());
			int x = a2b[v[k]];
			vv.push_back(x);
		}

		if(b == false) continue;

		for(int i = 0; i < vv.size(); i++) vv[i]--;
		hs.add_node_list(vv, c);
	}
	return 0;
}
