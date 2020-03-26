#include "essential.h"
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

bridger::bridger(splice_graph &g, const vector<PRC> &v)
	: gr(g), vpr(v)
{
	dp_solution_size = 10;
	dp_stack_size = 5;
	length_median = 300;
	length_low = 100;
	length_high = 500;
}

int bridger::resolve()
{
	build_piers();
	nominate();
	vote();
	return 0;
}

int bridger::build_piers()
{
	piers.clear();
	int n = gr.num_vertices();
	set<PI> ss;
	for(int k = 0; k < vpr.size(); k++)
	{
		if(vpr[k].first.vv.size() <= 0) continue;
		if(vpr[k].second.vv.size() <= 0) continue;
		int s = vpr[k].first.vv.back();
		int t = vpr[k].second.vv.front();
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

int bridger::build_piers_index()
{
	for(int k = 0; k < piers.size(); k++)
	{
		PI p(piers[k].bs, piers[k].bt);
		pindex.insert(pair<PI, int>(p, k));
	}
	return 0;
}

int bridger::vote()
{
	build_piers_index();
	opt.resize(vpr.size());
	for(int i = 0; i < vpr.size(); i++)
	{
		vote(vpr[i], opt[i]);
	}
	return 0;
}

int bridger::vote(const PRC &prc, phase &bbp)
{
	const rcluster &r1 = prc.first;
	const rcluster &r2 = prc.second;

	if(r1.vv.size() == 0) return 0;
	if(r2.vv.size() == 0) return 0;

	bbp.type = -1;
	int ss = r1.vv.back();
	int tt = r2.vv.front();

	// construct candidate bridging paths
	int type = 0;
	vector< vector<int> > pn;
	vector<int> ps;
	if(ss >= tt)
	{
		vector<int> v;
		bool b = merge_phasing_paths(r1.vv, r2.vv, v);
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
			vector<int> px = r1.vv;
			if(pb[e].v.size() >= 2) px.insert(px.end(), pb[e].v.begin() + 1, pb[e].v.end() - 1);
			px.insert(px.end(), r2.vv.begin(), r2.vv.end());
			pn.push_back(px);
			ps.push_back(pb[e].score);
		}
	}

	if(pn.size() == 0) return 0;

	vector<int> votes(pn.size(), 0);

	assert(r1.vl.size() == r2.vl.size());
	assert(r1.vr.size() == r2.vr.size());
	vector<int> bulls(r1.vl.size(), 0);

	vector<int32_t> upper_length;
	for(int e = 0; e < pn.size(); e++)
	{
		int32_t length = gr.get_total_length_of_vertices(pn[e]);
		upper_length.push_back(length);
	}

	for(int j = 0; j < r1.vl.size(); j++)
	{
		for(int e = 0; e < pn.size(); e++)
		{
			int32_t length = upper_length[e] - r1.vl[j] - r2.vr[j];
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

	if(votes[be] <= 0) return 0;
	if(voted <= 0) return 0;

	double voting_ratio = 100.0 * voted / r1.vl.size();
	double best_ratio = 100.0 * votes[be] / voted;

	// TODO parameters
	if(pn.size() >= 2 && ps[be] <= 1) return 0;
	if(voting_ratio <= 0.49) return 0;
	if(best_ratio < 0.8 && be != 0) return 0;

	bbp.type = type;
	bbp.score = ps[be];
	bbp.count = votes[be];
	bbp.v = pn[be];

	// TODO TODO
	// collect unbridged individuals
	/*
	for(int j = 0; j < r1.vl.size(); j++)
	{
		if(bulls[j] != be) continue;
		bridged[fr.h1] = true;
		bridged[fr.h2] = true;
	}
	*/

	/*
	   printf("fcluster %d: %lu fragments, %d voted, %lu candidates, best = %d, voted-best = %d, score = %d, type = %d, v = ( ", 
	   i, fc.frset.size(), voted, pn.size(), be, votes[be], ps[be], fc.bbp.type);
	   printv(fc.bbp.v);
	   printf("), v1 = ( ");
	   printv(r1.vv);
	   printf("), v2 = ( ");
	   printv(r2.vv);
	   printf(")\n");
	 */
	return 0;
}

int bridger::collect_unbridged_clusters(vector<PRC> &v)
{
	v.clear();
	for(int i = 0; i < opt.size(); i++)
	{
		if(opt[i].type >= 0) continue;
		v.push_back(vpr[i]);
	}
	return 0;
}

int bridger::build_hyper_set(hyper_set &hs)
{
	assert(opt.size() == vpr.size());
	for(int i = 0; i < opt.size(); i++)
	{
		if(opt[i].type >= 0)
		{
			vector<int> v = opt[i].v;
			for(int k = 0; k < v.size(); k++) v[k]--;
			hs.add_node_list(v, vpr[i].first.vl.size());
		}
		else
		{
			vector<int> v1 = vpr[i].first.vv;
			vector<int> v2 = vpr[i].second.vv;
			for(int k = 0; k < v1.size(); k++) v1[k]--;
			for(int k = 0; k < v2.size(); k++) v2[k]--;
			hs.add_node_list(v1, vpr[i].first.vl.size());
			hs.add_node_list(v2, vpr[i].second.vl.size());
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
