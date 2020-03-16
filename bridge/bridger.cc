/*
Part of Coral
(c) 2019 by Mingfu Shao, The Pennsylvania State University
See LICENSE for licensing.
*/

#include "bridger.h"
#include "util.h"
#include "config.h"

#include <algorithm>

int entry::print()
{
	printf("entry: length = %d, trace = (%d, %d), stack = (", length, trace1, trace2);
	printv(stack);
	printf(")\n");
	return 0;
}

bool entry_compare(const entry &x, const entry &y)
{
	for(int i = 0; i < x.stack.size() && i < y.stack.size(); i++)
	{
		if(x.stack[i] > y.stack[i]) return true;
		if(x.stack[i] < y.stack[i]) return false;
	}
	if(x.length < y.length) return true;
	else return false;
}

bridger::bridger(splice_graph &g, vector<hit> &h)
	: gr(g), hits(h)
{}

int bridger::resolve()
{
	build_fragments();
	build_piers();
	bridge();
	return 0;
}

int bridger::build_vertex_index()
{
	lindex.clear();
	rindex.clear();
	int n = gr.num_vertices() - 1;
	for(int i = 0; i <= n; i++)
	{
		PI32 p = regions[i].first;
		const vertex_info &v = gr.get_vertex_info(i);
		if(i != 0) lindex.insert(pair<int32_t, int>(v.lpos, i));
		if(i != n) rindex.insert(pair<int32_t, int>(v.rpos, i));
	}
	return 0;
}

int bridger::build_fragments()
{
	// TODO: parameters
	int32_t max_misalignment1 = 20;
	int32_t max_misalignment2 = 10;

	vector<bool> paired(hits.size(), false);

	fragments.clear();
	if(hits.size() == 0) return 0;

	int max_index = hits.size() + 1;
	if(max_index > 1000000) max_index = 1000000;

	vector< vector<int> > vv;
	vv.resize(max_index);

	// first build index
	for(int i = 0; i < hits.size(); i++)
	{
		hit &h = hits[i];
		if(h.isize >= 0) continue;
		//if(h.vlist.size() == 0) continue;

		// do not use hi; as long as qname, pos and isize are identical
		int k = (h.get_qhash() % max_index + h.pos % max_index + (0 - h.isize) % max_index) % max_index;
		/*
		SI si(h.qname, h.hi);
		MSI &m = vv[k];
		assert(m.find(si) == m.end());
		m.insert(PSI(si, i));
		*/
		vv[k].push_back(i);
	}

	for(int i = 0; i < hits.size(); i++)
	{
		hit &h = hits[i];
		if(paired[i] == true) continue;
		if(h.isize <= 0) continue;
		//if(h.vlist.size() == 0) continue;

		int k = (h.get_qhash() % max_index + h.mpos % max_index + h.isize % max_index) % max_index;

		/*
		h.print();
		for(int j = 0; j < vv[k].size(); j++)
		{
			hit &z = hits[vv[k][j]];
			printf(" ");
			z.print();
		}
		*/

		int x = -1;
		for(int j = 0; j < vv[k].size(); j++)
		{
			int u = vv[k][j];
			hit &z = hits[u];
			//if(z.hi != h.hi) continue;
			if(paired[u] == true) continue;
			if(z.pos != h.mpos) continue;
			if(z.isize + h.isize != 0) continue;
			//if(z.qhash != h.qhash) continue;
			if(z.qname != h.qname) continue;
			x = vv[k][j];
			break;
		}

		/*
		SI si(h.qname, h.hi);
		MSI::iterator it = vv[k].find(si);
		if(it == vv[k].end()) continue;
		int x = it->second;
		*/

		//printf("HIT: i = %d, x = %d, hits[i].vlist = %lu | ", i, x, hits[i].vlist.size(), hits[i].qname.c_str()); hits[i].print();

		if(x == -1) continue;
		//if(hits[x].vlist.size() == 0) continue;

		fragment fr(&hits[i], &hits[x]);
		fr.lpos = h.pos;
		fr.rpos = hits[x].rpos;

		/*
		vector<int> v1 = decode_vlist(hits[i].vlist);
		vector<int> v2 = decode_vlist(hits[x].vlist);
		fr.k1l = fr.h1->pos - regions[v1.front()].lpos;
		fr.k1r = regions[v1.back()].rpos - fr.h1->rpos;
		fr.k2l = fr.h2->pos - regions[v2.front()].lpos;
		fr.k2r = regions[v2.back()].rpos - fr.h2->rpos;

		fr.b1 = true;
		if(v1.size() <= 1) 
		{
			fr.b1 = false;
		}
		else if(v1.size() >= 2 && v1[v1.size() - 2] == v1.back() - 1)
		{
			if(fr.h1->rpos - regions[v1.back()].lpos > max_misalignment1 + fr.h1->nm) fr.b1 = false;
		}
		else if(v1.size() >= 2 && v1[v1.size() - 2] != v1.back() - 1)
		{
			if(fr.h1->rpos - regions[v1.back()].lpos > max_misalignment2 + fr.h1->nm) fr.b1 = false;
		}

		fr.b2 = true;
		if(v2.size() <= 1)
		{
			fr.b2 = false;
		}
		else if(v2.size() >= 2 || v2[1] == v2.front() + 1)
		{
			if(regions[v2.front()].rpos - fr.h2->pos > max_misalignment1 + fr.h2->nm) fr.b2 = false;
		}
		else if(v2.size() >= 2 || v2[1] != v2.front() + 1)
		{
			if(regions[v2.front()].rpos - fr.h2->pos > max_misalignment2 + fr.h2->nm) fr.b2 = false;
		}
		*/

		fragments.push_back(fr);

		paired[i] = true;
		paired[x] = true;
	}

	//printf("total hits = %lu, total fragments = %lu\n", hits.size(), fragments.size());
	return 0;
}

int bridger::build_piers()
{
	piers.clear();
	map<PI, int> m;
	int n = gr.num_vertices();
	for(int k = 0; k < fragments.size(); k++)
	{
		fragment &fr = fragments[k];
		fr.pr = NULL;
		int32_t p = fr.h1->rpos;
		int32_t q = fr.h2->pos;
		int s = locate_vertex(p - 1, 0, n);
		int t = locate_vertex(q, 0, n);
		if(s < 0 || t < 0) continue;

		if(m.find(PI(s, t)) == m.end())
		{
			pier pr(s, t);
			m.insert(pair<PI, int>(PI(s, t), piers.size()));
			piers.push_back(pr);
			fr.pr = &(piers.back());
		}
		else
		{
			fr.pr = &(piers[m[PI(s, t)]]);
		}
	}
	return 0;
}

int bridger::bridge()
{
	if(piers.size() <= 0) return 0;

	sort(piers.begin(), piers.end());

	vector<int> bounds;
	bounds.push_back(0);
	for(int i = 1; i < piers.size(); i++)
	{
		if(piers[i].bs != piers[i - 1].bs)
		{
			bounds.push_back(i - 1);
			bounds.push_back(i - 0);
		}
	}
	bounds.push_back(piers.size() - 1);

	vector< vector<entry> > table;
	table.resize(gr.num_vertices());
	for(int k = 0; k < bounds.size() / 2; k++)
	{
		int b1 = bounds[k * 2 + 0];
		int b2 = bounds[k * 2 + 1];
		assert(piers[b1].bs == piers[b2].bs);
		int k1 = piers[b2].bs;
		int k2 = piers[b2].bt;

		dynamic_programming(k1, k2, table);

		for(int b = b1; b <= b2; b++)
		{
			int bt = piers[b].bt;
			vector< vector<int> > pb = trace_back(bt, table);

			for(int j = 0; j < pb.size(); j++)
			{
				path p;
				p.score = table[bt][j].stack.front();
				p.v = pb[j];
				piers[b].paths.push_back(p);
			}
		}
	}
	return 0;
}

int bridger::dynamic_programming(int k1, int k2, vector< vector<entry> > &table)
{
	int n = gr.num_vertices();
	assert(k1 >= 0 && k1 < n);
	assert(k2 >= 0 && k2 < n);

	table.clear();
	table.resize(n);

	table[k1].resize(1);
	table[k1][0].stack.assign(dp_stack_size, 999999);
	table[k1][0].length = gr.get_vertex_info(k1).rpos - gr.get_vertex_info(k1).lpos;
	table[k1][0].trace1 = -1;
	table[k1][0].trace2 = -1;

	for(int k = k1 + 1; k <= k2; k++)
	{
		vector<entry> v;
		int32_t len = gr.get_vertex_info(k).rpos - gr.get_vertex_info(k).lpos;
		PEEI pi = gr.in_edges(k);
		for(edge_iterator it = pi.first; it != pi.second; it++)
		{
			edge_descriptor e = (*it);
			int j = e->source();
			int w = (int)(gr.get_edge_weight(e));
			if(j < k1) continue;
			if(table[j].size() == 0) continue;

			for(int i = 0; i < table[j].size(); i++)
			{
				entry e;
				e.stack = update_stack(table[j][i].stack, w);
				e.length = table[j][i].length + len;
				e.trace1 = j;
				e.trace2 = i;
				v.push_back(e);
			}
		}

		sort(v.begin(), v.end(), entry_compare);
		if(v.size() > dp_solution_size) v.resize(dp_solution_size);
		table[k] = v;
	}
	return 0;
}

int bridger::locate_vertex(int32_t p, int a, int b)
{
	if(a >= b) return -1;
	int m = (a + b) / 2;
	assert(m >= 0 && m < gr.num_vertices());
	const vertex_info &v = gr.get_vertex_info(m);
	if(p >= v.lpos && p < v.rpos) return m;
	if(p < v.lpos) return locate_vertex(p, a, m - 1);
	return locate_vertex(p, m + 1, b);
}

vector<int> bridger::update_stack(const vector<int> &v, int s)
{
	vector<int> stack(v.size(), 0);
	for(int i = 0, j = 0; i < v.size() && j < v.size(); i++, j++)
	{
		if(i == j && v[i] > s)
		{
			stack[j] = s;
			j++;
		}
		stack[j] = v[i];
	}
	return stack;
}

vector< vector<int> > bridger::trace_back(int k, const vector< vector<entry> > &table)
{
	vector< vector<int> > vv;
	for(int i = 0; i < table[k].size(); i++)
	{
		vector<int> v;
		int p = k;
		int q = i;
		while(true)
		{
			v.push_back(p);
			const entry &e = table[p][q];
			p = e.trace1;
			q = e.trace2;
			if(p < 0) break;
		}
		reverse(v);
		vv.push_back(v);
	}
	return vv;
}

