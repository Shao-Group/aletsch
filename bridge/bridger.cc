/*
Part of Coral
(c) 2019 by Mingfu Shao, The Pennsylvania State University
See LICENSE for licensing.
*/

#include "hyper_graph.h"
#include "bridger.h"
#include "util.h"

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
{
	dp_solution_size = 10;
	dp_stack_size = 5;
	length_median = 300;
	length_low = 100;
	length_high = 500;
}

int bridger::resolve()
{
	//printf("start bridger: %lu hits, splice graph with %lu vertices and %lu edges\n", hits.size(), gr.num_vertices(), gr.num_edges());
	//gr.print();

	build_vertex_index();

	build_fragments();
	//printf("built %lu fragments\n", fragments.size());

	build_fclusters();
	//printf("built %lu fclusters\n", fclusters.size());

	build_piers();
	//printf("built %lu piers\n", piers.size());

	nominate();
	//printf("finish bridging\n");

	vote();
	//printf("finish voting\n\n");

	build_hyper_set();
	//printf("finish building hyper-set\n\n");

	return 0;
}

int bridger::build_vertex_index()
{
	lindex.clear();
	rindex.clear();
	int n = gr.num_vertices() - 1;
	for(int i = 0; i <= n; i++)
	{
		const vertex_info &v = gr.get_vertex_info(i);
		if(i != 0) lindex.insert(pair<int32_t, int>(v.lpos, i));
		if(i != n) rindex.insert(pair<int32_t, int>(v.rpos, i));
	}
	return 0;
}

bool bridger::align_hit(const hit &h, vector<int> &vv)
{
	vv.clear();
	vector<int64_t> v;
	h.get_aligned_intervals(v);
	if(v.size() == 0) return false;

	vector<PI> sp;
	sp.resize(v.size());

	int32_t p1 = high32(v.front());
	int32_t p2 = low32(v.back());

	sp[0].first = locate_vertex(p1, 0, gr.num_vertices());
	//if(sp[0].first < 0) printf("start boundary fail\n");
	if(sp[0].first < 0) return false;

	for(int k = 1; k < v.size(); k++)
	{
		p1 = high32(v[k]);
		map<int32_t, int>::const_iterator it = lindex.find(p1);
		//if(it == lindex.end()) printf("lindex fail\n");
		if(it == lindex.end()) return false;
		sp[k].first = it->second;
	}

	sp[sp.size() - 1].second = locate_vertex(p2 - 1, 0, gr.num_vertices());
	//if(sp[sp.size() - 1].second < 0) printf("end boundary fail\n");
	if(sp[sp.size() - 1].second < 0) return false;

	for(int k = 0; k < v.size() - 1; k++)
	{
		p2 = low32(v[k]);
		map<int32_t, int>::const_iterator it = rindex.find(p2);
		//if(it == rindex.end()) printf("rindex fail\n");
		if(it == rindex.end()) return false;
		sp[k].second = it->second; 
	}

	for(int k = 0; k < sp.size(); k++)
	{
		assert(sp[k].first <= sp[k].second);
		if(k > 0) assert(sp[k - 1].second < sp[k].first);
		for(int j = sp[k].first; j <= sp[k].second; j++) vv.push_back(j);
	}
	return true;
}

int bridger::build_fragments()
{
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

		// do not use hi; as long as qname, pos and isize are identical
		int k = (h.get_qhash() % max_index + h.pos % max_index + (0 - h.isize) % max_index) % max_index;
		vv[k].push_back(i);
	}

	for(int i = 0; i < hits.size(); i++)
	{
		hit &h = hits[i];
		if(paired[i] == true) continue;
		if(h.isize <= 0) continue;

		int k = (h.get_qhash() % max_index + h.mpos % max_index + h.isize % max_index) % max_index;
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

		if(x == -1) continue;

		fragment fr(i, x);
		fr.lpos = h.pos;
		fr.rpos = hits[x].rpos;

		fragments.push_back(fr);

		paired[i] = true;
		paired[x] = true;
	}

	//printf("total hits = %lu, total fragments = %lu\n", hits.size(), fragments.size());
	return 0;
}

int bridger::build_fclusters()
{
	typedef pair< vector<int>, vector<int> > PVV;
	map<PVV, int> findex;			// index for fclusters
	fclusters.clear();

	// TODO: parameters
	int32_t max_misalignment1 = 20;
	int32_t max_misalignment2 = 10;

	int aligned = 0;
	for(int i = 0; i < fragments.size(); i++)
	{
		fragment &fr = fragments[i];

		vector<int> v1;
		vector<int> v2;
		bool b1 = align_hit(hits[fr.h1], v1);
		bool b2 = align_hit(hits[fr.h2], v2);
		if(b1 == false || b2 == false) continue;

		aligned++;

		PVV pvv(v1, v2);
		if(findex.find(pvv) == findex.end())
		{
			fcluster fc;
			fc.v1 = v1;
			fc.v2 = v2;
			fc.frset.push_back(i);
			findex.insert(pair<PVV, int>(pvv, fclusters.size()));
			fclusters.push_back(fc);
		}
		else
		{
			int k = findex[pvv];
			fcluster &fc = fclusters[k];
			assert(fc.v1 == v1);
			assert(fc.v2 == v2);
			fc.frset.push_back(i);
		}

		// setup fragment
		fr.k1l = hits[fr.h1].pos - gr.get_vertex_info(v1.front()).lpos;
		fr.k1r = gr.get_vertex_info(v1.back()).rpos - hits[fr.h1].rpos;
		fr.k2l = hits[fr.h2].pos - gr.get_vertex_info(v2.front()).lpos;
		fr.k2r = gr.get_vertex_info(v2.back()).rpos - hits[fr.h2].rpos;

		fr.b1 = true;
		if(v1.size() <= 1) 
		{
			fr.b1 = false;
		}
		else if(v1.size() >= 2 && v1[v1.size() - 2] == v1.back() - 1)
		{
			if(hits[fr.h1].rpos - gr.get_vertex_info(v1.back()).lpos > max_misalignment1 + hits[fr.h1].nm) fr.b1 = false;
		}
		else if(v1.size() >= 2 && v1[v1.size() - 2] != v1.back() - 1)
		{
			if(hits[fr.h1].rpos - gr.get_vertex_info(v1.back()).lpos > max_misalignment2 + hits[fr.h1].nm) fr.b1 = false;
		}

		fr.b2 = true;
		if(v2.size() <= 1)
		{
			fr.b2 = false;
		}
		else if(v2.size() >= 2 || v2[1] == v2.front() + 1)
		{
			if(gr.get_vertex_info(v2.front()).rpos - hits[fr.h2].pos > max_misalignment1 + hits[fr.h2].nm) fr.b2 = false;
		}
		else if(v2.size() >= 2 || v2[1] != v2.front() + 1)
		{
			if(gr.get_vertex_info(v2.front()).rpos - hits[fr.h2].pos > max_misalignment2 + hits[fr.h2].nm) fr.b2 = false;
		}
	}

	//printf(" total %d / %lu aligned fragments\n", aligned, fragments.size());
	return 0;
}

int bridger::build_piers()
{
	piers.clear();
	int n = gr.num_vertices();
	set<PI> ss;
	for(int k = 0; k < fclusters.size(); k++)
	{
		fcluster &fc = fclusters[k];
		if(fc.v1.size() <= 0) continue;
		if(fc.v2.size() <= 0) continue;
		int s = fc.v1.back();
		int t = fc.v2.front();
		if(s >= t) continue;
		ss.insert(PI(s, t));
	}

	for(set<PI>::iterator it = ss.begin(); it != ss.end(); it++)
	{
		pier pr(it->first, it->second);
		piers.push_back(pr);
	}
	return 0;
}

int bridger::nominate()
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
				phase p;
				p.score = table[bt][j].stack.front();
				p.stack = table[bt][j].stack;
				p.v = pb[j];
				piers[b].phases.push_back(p);
			}

			sort(piers[b].phases.begin(), piers[b].phases.end(), compare_phase_stack);

			/*
			printf("phases for pier %d: (%d, %d)\n", b, piers[b].bs, piers[b].bt);
			for(int i = 0; i < piers[b].phases.size(); i++)
			{
				piers[b].phases[i].print(i);
			}
			*/
		}
	}
	return 0;
}

int bridger::vote()
{
	bridged.resize(hits.size(), false);

	// build index for piers
	map<PI, int> pindex;
	for(int k = 0; k < piers.size(); k++)
	{
		PI p(piers[k].bs, piers[k].bt);
		pindex.insert(pair<PI, int>(p, k));
	}

	for(int i = 0; i < fclusters.size(); i++)
	{
		fcluster &fc = fclusters[i];
		if(fc.v1.size() == 0) continue;
		if(fc.v2.size() == 0) continue;

		fc.bbp.type = -1;
		int ss = fc.v1.back();
		int tt = fc.v2.front();

		// construct candidate bridging paths
		int type = 0;
		vector< vector<int> > pn;
		vector<int> ps;
		if(ss >= tt)
		{
			vector<int> v;
			bool b = merge_phasing_paths(fc.v1, fc.v2, v);
			if(b == true) 
			{
				type = 1;
				pn.push_back(v);
				ps.push_back(10);
			}
		}
		else if(pindex.find(PI(ss, tt)) != pindex.end())
		{
			type = 2;
			int k = pindex[PI(ss, tt)];
			vector<phase> &pb = piers[k].phases;
			for(int e = 0; e < pb.size(); e++)
			{
				vector<int> px = fc.v1;
				if(pb[e].v.size() >= 2) px.insert(px.end(), pb[e].v.begin() + 1, pb[e].v.end() - 1);
				px.insert(px.end(), fc.v2.begin(), fc.v2.end());
				pn.push_back(px);
				ps.push_back(pb[e].score);
			}
		}

		if(pn.size() == 0) continue;

		vector<int> votes(pn.size(), 0);
		vector<int> bulls(fc.frset.size(), -1);
		for(int j = 0; j < fc.frset.size(); j++)
		{
			fragment &fr = fragments[fc.frset[j]];
			for(int e = 0; e < pn.size(); e++)
			{
				int32_t length = compute_aligned_length(fr, pn[e]);
				if(length < length_low) continue;
				if(length > length_high) continue;
				votes[e]++;
				bulls[j] = e;
				break;
			}
		}

		int be = 0;
		int voted = votes[0];
		for(int j = 1; j < votes.size(); j++)
		{
			voted += votes[j];
			if(votes[j] > votes[be]) be = j;
		}

		if(votes[be] <= 0) continue;
		if(voted <= 0) continue;

		double voting_ratio = 100.0 * voted / fc.frset.size();
		double best_ratio = 100.0 * votes[be] / voted;

		// TODO parameters
		/*
		if(voting_ratio <= 0.49) continue;
		if(best_ratio < 0.8 && be != best_phase) continue;
		*/

		for(int j = 0; j < fc.frset.size(); j++)
		{
			fragment &fr = fragments[fc.frset[j]];
			if(bulls[j] != be) continue;
			bridged[fr.h1] = true;
			bridged[fr.h2] = true;
		}

		fc.bbp.type = type;
		fc.bbp.score = ps[be];
		fc.bbp.count = votes[be];
		fc.bbp.v = pn[be];

		/*
		printf("fcluster %d: %lu fragments, %d voted, %lu candidates, best = %d, voted-best = %d, score = %d, type = %d, v = ( ", 
				i, fc.frset.size(), voted, pn.size(), be, votes[be], ps[be], fc.bbp.type);
		printv(fc.bbp.v);
		printf("), v1 = ( ");
		printv(fc.v1);
		printf("), v2 = ( ");
		printv(fc.v2);
		printf(")\n");
		*/
	}
	return 0;
}

int bridger::build_hyper_set()
{
	hs.clear();

	int c1 = 0;
	for(int i = 0; i < fclusters.size(); i++)
	{
		fcluster &fc = fclusters[i];
		if(fc.bbp.type < 0) continue;
		int c = fc.bbp.count;
		vector<int> v = fc.bbp.v;
		for(int k = 0; k < v.size(); k++) v[k]--;
		hs.add_node_list(v, c);
		c1 += c * 2;
	}

	int c2 = 0;
	for(int i = 0; i < hits.size(); i++)
	{
		if(bridged[i] == true) continue;
		vector<int> v;
		bool b = align_hit(hits[i], v);
		if(b == false) continue;
		for(int k = 0; k < v.size(); k++) v[k]--;
		hs.add_node_list(v, 1);
		c2++;
	}

	//printf("total %d hits are paired, %d hits are isolated\n", c1, c2);
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
	if(p < v.lpos) return locate_vertex(p, a, m);
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

int32_t bridger::compute_aligned_length(fragment &fr, const vector<int>& v)
{
	if(v.size() == 0) return -1;
	int32_t flen = 0;
	for(int i = 0; i < v.size(); i++)
	{
		const vertex_info &x = gr.get_vertex_info(v[i]);
		flen += x.rpos - x.lpos;
	}
	return flen - fr.k1l - fr.k2r;
}
