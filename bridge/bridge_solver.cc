#include "essential.h"
#include "bridge_solver.h"
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

bridge_solver::bridge_solver(splice_graph &g, const vector<pereads_cluster> &v)
	: gr(g), vc(v)
{
	dp_solution_size = 10;
	dp_stack_size = 5;
	length_low = 100;
	length_high = 500;
	if(gr.lindex.size() == 0) 
	{
		assert(gr.rindex.size() == 0);
		gr.build_vertex_index();
	}
	build_bridging_vertices();
	build_piers();
	nominate();
	vote();
}

int bridge_solver::build_bridging_vertices()
{
	vpairs.clear();
	for(int i = 0; i < vc.size(); i++)
	{
		const pereads_cluster &pc = vc[i];
		int v1 = gr.locate_vertex(pc.bounds[1] - 1);
		int v2 = gr.locate_vertex(pc.bounds[1] - 0);
		vpairs.push_back(PI(v1, v2));
	}
	return 0;
}

int bridge_solver::build_piers()
{
	piers.clear();
	set<PI> ss;
	assert(vc.size() == vpairs.size());
	for(int k = 0; k < vc.size(); k++)
	{
		PI &p = vpairs[k];
		if(p.first < 0) continue;
		if(p.second < 0) continue;
		if(p.first >= p.second) continue;
		if(ss.find(p) != ss.end()) continue;
		ss.insert(p);
		pier pr(p.first, p.second);
		piers.push_back(pr);
	}
	return 0;
}

int bridge_solver::build_piers_index()
{
	pindex.clear();
	for(int k = 0; k < piers.size(); k++)
	{
		PI p(piers[k].bs, piers[k].bt);
		pindex.insert(pair<PI, int>(p, k));
	}
	return 0;
}

int bridge_solver::nominate()
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
				bridge_path p;
				p.score = table[bt][j].stack.front();
				p.stack = table[bt][j].stack;
				p.v = pb[j];
				build_intron_coordinates_from_path(gr, p.v, p.chain);
				piers[b].bridges.push_back(p);
			}

			sort(piers[b].bridges.begin(), piers[b].bridges.end(), compare_bridge_path_stack);

			/*
			printf("bridges for pier %d: (%d, %d)\n", b, piers[b].bs, piers[b].bt);
			for(int i = 0; i < piers[b].bridges.size(); i++)
			{
				piers[b].bridges[i].print(i);
			}
			*/
		}
	}
	return 0;
}

int bridge_solver::vote()
{
	build_piers_index();
	opt.resize(vc.size());
	for(int i = 0; i < vc.size(); i++)
	{
		vote(i, opt[i]);
	}
	return 0;
}

int bridge_solver::vote(int r, bridge_path &bbp)
{
	bbp.type = -1;
	int ss = vpairs[r].first;
	int tt = vpairs[r].second;
	if(ss < 0 || tt < 0) return 0;

	const pereads_cluster &pc = vc[r];

	// construct candidate bridging paths
	int type = 0;
	vector< vector<int32_t> > chains;
	vector<int> scores;
	if(ss >= tt)
	{
		vector<int32_t> v;
		if(pc.chain1.size() == 0 || pc.chain2.size() == 0 || pc.chain1.back() < pc.chain2.front())
		{
			type = 1;
			chains.push_back(v);
			scores.push_back(10);
		}
	}
	else if(pindex.find(PI(ss, tt)) != pindex.end())
	{
		type = 2;
		assert(pindex.find(PI(ss, tt)) != pindex.end());
		int k = pindex[PI(ss, tt)];
		vector<bridge_path> &pb = piers[k].bridges;
		for(int e = 0; e < pb.size(); e++)
		{
			chains.push_back(pb[e].chain);
			scores.push_back(pb[e].score);
		}
	}

	if(chains.size() == 0) return 0;

	int be = -1;
	int32_t length1 = get_total_length_of_introns(pc.chain1);
	int32_t length2 = get_total_length_of_introns(pc.chain2);
	for(int e = 0; e < chains.size(); e++)
	{
		if(pc.chain1.size() > 0 && chains[e].size() > 0) assert(pc.chain1.back() < chains[e].front());
		if(pc.chain2.size() > 0 && chains[e].size() > 0) assert(pc.chain2.front() > chains[e].back());

		int32_t length3 = get_total_length_of_introns(chains[e]);
		int32_t length = pc.bounds[3] - pc.bounds[0] - length1 - length2 - length3;
		if(length < length_low) continue;
		if(length > length_high) continue;
		be = e;
		break;
	}

	if(be < 0) return 0;

	bbp.type = type;
	bbp.score = scores[be];
	bbp.chain = chains[be];

	return 0;
}

int bridge_solver::collect_unbridged_clusters(vector<pereads_cluster> &v)
{
	v.clear();
	for(int i = 0; i < opt.size(); i++)
	{
		if(opt[i].type >= 0) continue;
		v.push_back(vc[i]);
	}
	return 0;
}

int bridge_solver::build_phase_set(phase_set &ps)
{
	assert(opt.size() == vc.size());
	for(int i = 0; i < vc.size(); i++)
	{
		add_phases_from_bridge_path(vc[i], opt[i], ps);
	}
	return 0;
}

int add_phases_from_bridge_path(const pereads_cluster &pc, const bridge_path &bbp, phase_set &ps)
{
	int32_t p0 = pc.extend[0];
	int32_t p3 = pc.extend[3];

	if(bbp.type >= 0)
	{
		vector<int32_t> v;
		v.push_back(p0);
		v.insert(v.end(), pc.chain1.begin(), pc.chain1.end());
		v.insert(v.end(), bbp.chain.begin(), bbp.chain.end());
		v.insert(v.end(), pc.chain2.begin(), pc.chain2.end());
		v.push_back(p3);
		ps.add(v, pc.count);
	}
	else
	{
		int32_t p1 = pc.extend[1];
		int32_t p2 = pc.extend[2];

		vector<int32_t> v1;
		v1.push_back(p0);
		v1.insert(v1.end(), pc.chain1.begin(), pc.chain1.end());
		v1.push_back(p1);
		ps.add(v1, pc.count);

		vector<int32_t> v2;
		v2.push_back(p2);
		v2.insert(v2.end(), pc.chain2.begin(), pc.chain2.end());
		v2.push_back(p3);

		ps.add(v2, pc.count);
	}
	return 0;
}

int bridge_solver::dynamic_programming(int k1, int k2, vector< vector<entry> > &table)
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

vector<int> bridge_solver::update_stack(const vector<int> &v, int s)
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

vector< vector<int> > bridge_solver::trace_back(int k, const vector< vector<entry> > &table)
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

int bridge_solver::print()
{
	assert(vc.size() == opt.size());
	int total_reads = 0;
	int bridged_reads = 0;
	int bridged_clusters = 0;
	for(int i = 0; i < vc.size(); i++)
	{
		total_reads += vc[i].count;
		if(opt[i].type < 0) continue;
		bridged_reads += vc[i].count;
		bridged_clusters++;
	}

	printf("bridge_solver: clusters %d / %lu, reads %d / %d, low = %d, high = %d\n", bridged_clusters, vc.size(), bridged_reads, total_reads, length_low, length_high);
	return 0;
}
