#include "util.h"
#include "essential.h"
#include "splice_graph.h"

#include <algorithm>

int transform_vertex_set_map(const set<int> &s, map<int, int> &m)
{
	m.clear();
	if(s.size() <= 0) return 0;

	vector<int> vv(s.begin(), s.end());
	sort(vv.begin(), vv.end());
	for(int i = 0; i < vv.size(); i++)
	{
		m.insert(pair<int, int>(vv[i], i + 1));
	}
	return 0;
}

int build_child_splice_graph(splice_graph &root, splice_graph &gr, map<int, int> &a2b)
{
	gr.clear();
	if(a2b.size() <= 0) return 0;

	vector<int> vv = get_keys(a2b);
	sort(vv.begin(), vv.end());

	gr.chrm = root.chrm;
	gr.strand = root.strand;

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

		//printf("there exists edge from 0 to %d in root\n", t);
		if(a2b.find(t) == a2b.end()) continue;
		int y = a2b[t];
		//printf("y = a2b[%d] = %d\n", t, y);

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
			assert(t == n || a2b.find(t) != a2b.end());
			int y = ((t == n) ? gr.num_vertices() - 1 : a2b[t]);

			edge_descriptor e = gr.add_edge(x, y);
			gr.set_edge_weight(e, root.get_edge_weight(*it1));
			gr.set_edge_info(e, root.get_edge_info(*it1));
		}
	}
	return 0;
}

int build_child_hyper_set(hyper_set &hyper, hyper_set &hs, map<int, int> &a2b)
{
	hs.clear();
	if(a2b.size() <= 0) return 0;

	for(MVII::const_iterator it = hyper.nodes.begin(); it != hyper.nodes.end(); it++)
	{
		vector<int> v = it->first;
		int c = it->second;

		if(v.size() <= 0) continue;
		if(a2b.find(v.front()) == a2b.end()) continue;

		// TODO: test
		vector<int> z = get_keys(a2b);
		printf("v = ( "); printv(v); printf(") \n");
		printf("z = ( "); printv(z); printf(") \n");

		vector<int> vv = project_vector(v, a2b);
		printf("V = ( "); printv(vv); printf(") \n");

		assert(vv.size() == v.size());

		for(int i = 0; i < vv.size(); i++) vv[i]--;
		hs.add_node_list(vv, c);
	}
	return 0;
}

vector<int> project_vector(const vector<int> &v, const map<int, int> &a2b)
{
	vector<int> vv;
	for(int k = 0; k < v.size(); k++)
	{
		map<int, int>::const_iterator it = a2b.find(v[k]);
		if(it == a2b.end()) break;
		vv.push_back(it->second);
	}
	return vv;
}

int build_path_coordinates(splice_graph &gr, const vector<int> &v, vector<int32_t> &vv)
{
	vv.clear();
	if(v.size() <= 0) return 0;

	int n = gr.num_vertices() - 1;
	int32_t pre = -99999;

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

	return 0;
}

bool check_continue_vertices(splice_graph &gr, int x, int y)
{
	if(x >= y) return true;
	for(int i = x; i < y; i++)
	{
		PEB p = gr.edge(i, i + 1);
		if(p.second == false) return false;
		if(gr.get_vertex_info(i).rpos != gr.get_vertex_info(i + 1).lpos) return false;
	}
	return true;
}

bool check_valid_path(splice_graph &gr, const vector<int> &vv)
{
	int n = gr.num_vertices() - 1;
	for(int k = 0; k < vv.size() - 1; k++)
	{
		if(vv[k + 0] < 0 && vv[k + 0] > n) return false;
		if(vv[k + 1] < 0 && vv[k + 1] > n) return false;
		PEB p = gr.edge(vv[k], vv[k + 1]);
		if(p.second == false) return false;
	}
	return true;
}
