/*
Part of Coral
(c) 2019 by Mingfu Shao, The Pennsylvania State University.
See LICENSE for licensing.
*/

#include "hyper_graph.h"
#include "config.h"
#include <algorithm>
#include <cstdio>

hyper_graph::hyper_graph(const MVII &m)
{
	nodes.clear();
	vcnts.clear();
	for(MVII::const_iterator it = m.begin(); it != m.end(); it++)
	{
		const vector<int> &v = it->first;
		int c = it->second;
		nodes.push_back(v);
		vcnts.push_back(c);
	}
}

int hyper_graph::build_overlap_index()
{
	adjsetx.resize(nodes.size());
	adjsety.resize(nodes.size());
	for(int k = 0; k < nodes.size(); k++) adjsetx[k].clear();
	for(int k = 0; k < nodes.size(); k++) adjsety[k].clear();

	for(int i = 0; i < nodes.size(); i++)
	{
		const vector<int> &vx = nodes[i];
		for(int j = i + 1; j < nodes.size(); j++)
		{
			const vector<int> &vy = nodes[j];

			PI p = get_overlap(vx, vy);
			if(p.first < 0 || p.second < 0) continue;
			assert(p.first >= 0 && p.second >= 0);

			adjsetx[i].insert(pair<int, int>(j, p.second));
			adjsety[j].insert(pair<int, int>(i, p.first));
		}
	}
	return 0;
}

PI hyper_graph::get_overlap(const vector<int> &vx, const vector<int> &vy) 
{
	if(vx.size() == 0) return PI(-1, -1);
	if(vy.size() == 0) return PI(-1, -1);

	// TODO TODO
	//if(vy.back() <= vx.back()) return -1;

	vector<int>::const_iterator ix = lower_bound(vx.begin(), vx.end(), vy[0]);
	if(ix == vx.end()) return PI(-1, -1);

	vector<int>::const_iterator iy = lower_bound(vy.begin(), vy.end(), vx.back());
	if(iy == vy.end()) return PI(-1, -1);

	int kx = ix - vx.begin();
	int ky = iy - vy.begin();

	if(vx.size() - kx != ky + 1) return PI(-1, -1);

	// TODO parameter, setup to 0
	double min_overlap_ratio = 0.6;
	if(vx.size() <= vy.size() && vx.size() - kx < min_overlap_ratio * (vx.size() - 1)) return PI(-1, -1);
	if(vy.size() <= vx.size() && ky + 1 < min_overlap_ratio * (vy.size() - 1)) return PI(-1, -1);

	bool b = identical(vx, kx, vx.size() - 1, vy, 0, ky);
	if(b == false) return PI(-1, -1);
	else return PI(kx, ky);
}

int hyper_graph::align_paths(const MVII &m)
{
	for(MVII::const_iterator it = m.begin(); it != m.end(); it++)
	{
		const vector<int> &v = it->first;
		int s = align_path(v);
		printf("align path, span = %d: list = ( ", s);
		printv(v);
		printf(")\n");
	}
	return 0;
}

int hyper_graph::align_path(const vector<int> &x)
{
	vector<int> open_node;
	vector<int> open_span;
	vector<int> open_node_index;
	vector<int> open_path_index;
	vector<bool> closed(nodes.size(), false);

	for(int i = 0; i < nodes.size(); i++)
	{
		vector<int>::const_iterator r = lower_bound(nodes[i].begin(), nodes[i].end(), x.front());
		if(r == nodes[i].end()) continue;
		if(*r != x.front()) continue;

		open_node.push_back(i);
		open_span.push_back(0);
		open_node_index.push_back(r - nodes[i].begin());
		open_path_index.push_back(0);
		closed[i] = true;
	}

	/* TODO
	map<int, int> m = index_path(x);
	for(int i = 0; i < adjsety.size(); i++)
	{
		if(closed[i] == true) continue;
		if(adjsety[i].size() >= 1) continue;
		int p = nodes[i].front();
		if(m.find(p) == m.end()) continue;
		int k = m[p];
		assert(k >= 1);
		assert(x[k] == p);

		open_node.push_back(i);
		open_span.push_back(1);
		open_node_index.push_back(0);
		open_path_index.push_back(k);
		closed[i] = true;
	}
	*/

	int q = 0;
	while(q < open_node.size())
	{
		int z = open_node[q];
		int s = open_span[q];
		int r = open_node_index[q];
		int k = open_path_index[q];

		q++;

		int d = x.size() - k <= nodes[z].size() - r ? x.size() - k : nodes[z].size() - r;
		bool b = identical(x, k, k + d - 1, nodes[z], r, r + d - 1);

		//for(int i = k; i <= k + d - 1; i++) printf("%d, ", x[i]); printf("\n");
		//for(int i = r; i <= r + d - 1; i++) printf("%d, ", nodes[z][i]); printf("\n");

		//printf(" check %d vertex: node = %d, node-index = %d, path-index = %d, d = %d: compatible = %d\n", q - 1, z, r, k, d, b ? 1 : 0);

		if(b == false) continue;
		if(d == x.size() - k) return s + 1;

		assert(d == nodes[z].size() - r);
		k += nodes[z].size() - r;

		// TODO
		//if(adjsetx[z].size() == 0) return s + 2;

		for(map<int, int>::iterator it = adjsetx[z].begin(); it != adjsetx[z].end(); it++)
		{
			int u = it->first;
			int j = it->second + 1;
			if(closed[u] == true) continue;
			if(j >= nodes[u].size()) continue;
			if(nodes[u][j] != x[k]) continue;

			open_node.push_back(u);
			open_span.push_back(s + 1);
			open_node_index.push_back(j);
			open_path_index.push_back(k);
			closed[u] = true;
		}
	}

	return -1;
}

map<int, int> hyper_graph::index_path(const vector<int> &x)
{
	map<int, int> m;
	for(int k = 0; k < x.size(); k++)
	{
		m.insert(pair<int, int>(x[k], k));
	}
	return m;
}

int hyper_graph::keep_compatible_nodes(splice_graph &gr)
{
	vector< vector<int> > zz;
	vector<int> uu;
	for(int i = 0; i < nodes.size(); i++)
	{
		int c = vcnts[i];
		const vector<int> &vv = nodes[i];
		if(vv.size() <= 1) continue;

		bool b = true;
		bool j = false;
		for(int k = 0; k < vv.size() - 1; k++)
		{
			PEB p = gr.edge(vv[k], vv[k + 1]);
			if(p.second == false) b = false;
			if(b == false) break;
			int32_t k1 = gr.get_vertex_info(vv[k + 0]).rpos;
			int32_t k2 = gr.get_vertex_info(vv[k + 1]).lpos;
			if(k1 + 1 < k2) j = true;
		}
		if(b == false) continue;
		//if(j == false) continue;
		zz.push_back(vv);
		uu.push_back(c);
	}
	nodes = zz;
	vcnts = uu;
	return 0;
}

int hyper_graph::keep_maximal_nodes()
{
	vector<bool> bb(nodes.size(), false);

	// build index with the first element
	map< int, vector<int> > m;
	for(int i = 0; i < nodes.size(); i++)
	{
		if(nodes[i].size() <= 0) continue;
		int s = nodes[i].front();

		if(m.find(s) == m.end())
		{
			vector<int> v;
			v.push_back(i);
			m.insert(pair<int, vector<int> >(s, v));
		}
		else
		{
			m[s].push_back(i);
		}
	}

	for(int i = 0; i < nodes.size(); i++)
	{
		if(bb[i] == true) continue;
		vector<int> &v = nodes[i];

		for(int a = 0; a < v.size(); a++)
		{
			int s = v[a];
			if(m.find(s) == m.end()) continue;

			vector<int> &x = m[s];
			for(int k = 0; k < x.size(); k++)
			{
				int j = x[k];
				vector<int> &u = nodes[j];

				// check whether i covers j
				assert(u.front() == s);
				if(j == i) continue;
				if(v.size() - a < u.size()) continue;
				bool b = identical(v, a, a + u.size() - 1, u, 0, u.size() - 1);
				if(b == true) bb[j] = true;
			}
		}
	}

	vector< vector<int> > zz;
	vector<int> uu;
	for(int i = 0; i < nodes.size(); i++)
	{
		if(bb[i] == true) continue;
		zz.push_back(nodes[i]);
		uu.push_back(vcnts[i]);
	}

	nodes = zz;
	vcnts = uu;

	return 0;
}

int hyper_graph::print_nodes()
{
	for(int i = 0; i < nodes.size(); i++)
	{
		printf("node %d, c = %d, list = ( ", i, vcnts[i]);
		printv(nodes[i]);
		printf(")\n");
	}
	return 0;
}

int hyper_graph::print_index()
{
	for(int i = 0; i < adjsetx.size(); i++)
	{
		for(map<int, int>::iterator it = adjsetx[i].begin(); it != adjsetx[i].end(); it++)
		{
			int j = it->first;
			int ky = it->second;
			int kx = adjsety[j][i];
			printf("overlap %d -> %d, index = (%d, %d)\n", i, j, kx, ky);
		}
	}
	return 0;
}

int compare_phasing_paths(const vector<int> &ref, const vector<int> &qry)
{
	if(ref.back() < qry.front()) return FALL_RIGHT;
	if(ref.front() > qry.back()) return FALL_LEFT;

	vector<int>::const_iterator r1 = lower_bound(ref.begin(), ref.end(), qry.front());
	vector<int>::const_iterator q1 = lower_bound(qry.begin(), qry.end(), ref.front());
	assert(r1 != ref.end());
	assert(q1 != qry.end());

	int kr1 = r1 - ref.begin();
	int kq1 = q1 - qry.begin();
	assert(kr1 == 0 || kq1 == 0);

	vector<int>::const_iterator q2 = lower_bound(qry.begin(), qry.end(), ref.back());
	vector<int>::const_iterator r2 = lower_bound(ref.begin(), ref.end(), qry.back());
	assert(r2 != ref.end() || q2 != qry.end());

	int kr2 = r2 - ref.begin();
	int kq2 = q2 - qry.begin();

	if(*q1 == ref.front() || *r1 == qry.front())
	{
		if(r2 != ref.end() && q2 != qry.end())
		{
			assert(ref.back() == qry.back());
			assert(kr2 == ref.size() - 1);
			assert(kq2 == qry.size() - 1);
			bool b = identical(ref, kr1, kr2, qry, kq1, kq2);
			if(b == false) return CONFLICTING;
			if(b == true && kr1 == 0 && kq1 == 0) return IDENTICAL;
			if(b == true && kr1 >= 1 && kq1 == 0) return CONTAINED;
			if(b == true && kr1 == 0 && kq1 >= 1) return CONTAINING;
			assert(false);
		}
		else if(r2 != ref.end() && q2 == qry.end())
		{
			bool b = identical(ref, kr1, kr2, qry, kq1, qry.size() - 1);
			if(b == false) return CONFLICTING;
			if(b == true && kq1 == 0) return CONTAINED;
			if(b == true && kq1 >= 1) return EXTEND_LEFT;
			assert(false);
		}
		else if(r2 == ref.end() && q2 != qry.end())
		{
			bool b = identical(ref, kr1, ref.size() - 1, qry, kq1, kq2);
			if(b == false) return CONFLICTING;
			if(b == true && kr1 == 0) return CONTAINING;
			if(b == true && kr1 >= 1) return EXTEND_RIGHT;
		}
	}
	else if(*r1 > qry.front() && r2 == r1 && *r2 > qry.back()) return NESTED;
	else if(*q1 > ref.front() && q2 == q1 && *q2 > ref.back()) return NESTING;
	return CONFLICTING;
}

bool identical(const vector<int> &x, int x1, int x2, const vector<int> &y, int y1, int y2)
{
	assert(x1 >= 0 && x1 < x.size());
	assert(x2 >= 0 && x2 < x.size());
	assert(y1 >= 0 && y1 < y.size());
	assert(y2 >= 0 && y2 < y.size());

	if(x[x1] != y[y1]) return false;
	if(x[x2] != y[y2]) return false;
	if(x2 - x1 != y2 - y1) return false;

	for(int kx = x1, ky = y1; kx <= x2 && ky <= y2; kx++, ky++)
	{
		if(x[kx] != y[ky]) return false;
	}

	return true;
}
