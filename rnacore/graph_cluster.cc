/*
Part of aletsch
(c) 2020 by Mingfu Shao, The Pennsylvania State University
See LICENSE for licensing.
*/

#include "graph_cluster.h"
#include "essential.h"
#include "util.h"

#include <algorithm>

graph_cluster::graph_cluster(splice_graph &g, bundle &d, int max_gap, bool b)
	: gr(g), bd(d), max_partition_gap(max_gap), store_hits(b)
{
	group_pereads();
} 

int graph_cluster::build_pereads_clusters(vector<pereads_cluster> &vc)
{
	for(int k = 0; k < groups.size(); k++)
	{
		build_pereads_clusters(k, vc);
	}
	return 0;
}

int graph_cluster::group_pereads()
{
	typedef pair< vector<int>, vector<int> > PVV;
	map<PVV, int> findex;

	extend.clear();
	groups.clear();
	for(int i = 0; i < bd.frgs.size(); i++)
	{
		// only group unbridged fragments
		if(bd.frgs[i][2] >= 1) continue;
		if(bd.frgs[i][2] <= -1) continue;

		bd.frgs[i][2] = -1;		// assume cannot be bridged

		int h1 = bd.frgs[i][0];
		int h2 = bd.frgs[i][1];

		assert(bd.hits[h1].hid >= 0);
		assert(bd.hits[h2].hid >= 0);

		if(bd.hits[h1].pos > bd.hits[h2].pos) continue;
		if(bd.hits[h1].rpos > bd.hits[h2].rpos) continue;

		vector<int> v1;
		vector<int> v2;
		vector<int32_t> chain1 = bd.hcst.get(h1).first;
		vector<int32_t> chain2 = bd.hcst.get(h2).first;

		bool b1 = align_hit_to_splice_graph(bd.hits[h1], chain1, gr, v1);
		bool b2 = align_hit_to_splice_graph(bd.hits[h2], chain2, gr, v2);

		if(b1 == false || b2 == false)  continue;
		if(v1.size() == 0 || v2.size() == 0) continue;

		bd.frgs[i][2] = 0;			// to be bridged

		PVV pvv(v1, v2);
		if(findex.find(pvv) == findex.end())
		{
			vector<int> v;
			v.push_back(i);
			findex.insert(pair<PVV, int>(pvv, groups.size()));
			int32_t p1 = gr.get_vertex_info(v1.front()).lpos;
			int32_t p2 = gr.get_vertex_info(v1.back()).rpos;
			int32_t p3 = gr.get_vertex_info(v2.front()).lpos;
			int32_t p4 = gr.get_vertex_info(v2.back()).rpos;
			extend.push_back(p1);
			extend.push_back(p2);
			extend.push_back(p3);
			extend.push_back(p4);
			groups.push_back(v);
		}
		else
		{
			int k = findex[pvv];
			groups[k].push_back(i);
		}
	}

	//for(int k = 0; k < groups.size(); k++) printf("group %d contains %lu frags\n", k, groups[k].size());

	return 0;
}

int graph_cluster::build_pereads_clusters(int g, vector<pereads_cluster> &vc)
{
	const vector<int> &fs = groups[g];
	vector<vector<int32_t>> vv;
	for(int i = 0; i < fs.size(); i++)
	{
		int h1 = bd.frgs[fs[i]][0];
		int h2 = bd.frgs[fs[i]][1];

		vector<int32_t> v;
		v.push_back(bd.hits[h1].pos);
		v.push_back(bd.hits[h1].rpos);
		v.push_back(bd.hits[h2].pos);
		v.push_back(bd.hits[h2].rpos);
		v.push_back(i);
		vv.push_back(v);
	}

	vector< vector<int> > zz = partition(vv, 0);

	for(int i = 0; i < zz.size(); i++)
	{
		if(zz[i].size() == 0) continue;

		int h1 = bd.frgs[fs[zz[i][0]]][0];
		int h2 = bd.frgs[fs[zz[i][0]]][1];
		assert(bd.hits[h1].rpos <= bd.hits[h2].rpos);
		assert(bd.hits[h1].pos <= bd.hits[h2].pos);

		pereads_cluster pc;
		pc.count = 0;
		pc.chain1 = bd.hcst.get(h1).first;
		pc.chain2 = bd.hcst.get(h2).first;

		vector<int32_t> bounds(4, 0);
		bounds[0] = bd.hits[h1].pos;
		bounds[1] = bd.hits[h1].rpos;
		bounds[2] = bd.hits[h2].pos;
		bounds[3] = bd.hits[h2].rpos;

		for(int k = 0; k < zz[i].size(); k++)
		{
			h1 = bd.frgs[fs[zz[i][k]]][0];
			h2 = bd.frgs[fs[zz[i][k]]][1];

			pc.bounds[0] += bd.hits[h1].pos  - bounds[0];
			pc.bounds[1] += bd.hits[h1].rpos - bounds[1];
			pc.bounds[2] += bd.hits[h2].pos  - bounds[2];
			pc.bounds[3] += bd.hits[h2].rpos - bounds[3];

			pc.frlist.push_back(fs[zz[i][k]]);
			pc.count++;
			
			if(store_hits == true)
			{
				pc.hits1.push_back(bd.hits[h1]);
				pc.hits2.push_back(bd.hits[h2]);
			}
		}

		if(pc.count <= 0) continue;

		pc.bounds[0] = pc.bounds[0] / pc.count + bounds[0];  
		pc.bounds[1] = pc.bounds[1] / pc.count + bounds[1];
		pc.bounds[2] = pc.bounds[2] / pc.count + bounds[2]; 
		pc.bounds[3] = pc.bounds[3] / pc.count + bounds[3];
		pc.extend[0] = extend[g * 4 + 0];
		pc.extend[1] = extend[g * 4 + 1];
		pc.extend[2] = extend[g * 4 + 2];
		pc.extend[3] = extend[g * 4 + 3];

		vc.push_back(pc);
	}

	return 0;
}

vector< vector<int> > graph_cluster::partition(vector< vector<int32_t> > &fs, int r)
{
	vector< vector<int> > vv;
	if(fs.size() == 0) return vv;

	if(r >= 4)
	{
		vector<int> v;
		for(int k = 0; k < fs.size(); k++) v.push_back(fs[k][4]);
		vv.push_back(v);
		return vv;
	}

	if(r == 0) sort(fs.begin(), fs.end(), compare_rank0);	
	if(r == 1) sort(fs.begin(), fs.end(), compare_rank1);	
	if(r == 2) sort(fs.begin(), fs.end(), compare_rank2);	
	if(r == 3) sort(fs.begin(), fs.end(), compare_rank3);	

	int pre = 0;
	for(int k = 1; k <= fs.size(); k++)
	{
		if(k < fs.size()) assert(fs[k][r] >= fs[k - 1][r]);
		if(k < fs.size() && fs[k][r] - fs[k - 1][r] <= max_partition_gap) continue;

		vector< vector<int32_t> > fs1;
		for(int i = pre; i < k; i++) fs1.push_back(fs[i]);

		vector< vector<int> > vv1 = partition(fs1, r + 1);
		vv.insert(vv.end(), vv1.begin(), vv1.end());

		pre = k;
	}
	return vv;
}

bool compare_rank0(const vector<int32_t> &x, const vector<int32_t> &y) { return x[0] < y[0]; }
bool compare_rank1(const vector<int32_t> &x, const vector<int32_t> &y) { return x[1] < y[1]; }
bool compare_rank2(const vector<int32_t> &x, const vector<int32_t> &y) { return x[2] < y[2]; }
bool compare_rank3(const vector<int32_t> &x, const vector<int32_t> &y) { return x[3] < y[3]; }
