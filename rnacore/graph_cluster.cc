/*
Part of meta-scallop
(c) 2020 by Mingfu Shao, The Pennsylvania State University
See LICENSE for licensing.
*/


#include "graph_cluster.h"
#include "essential.h"
#include "util.h"

#include <algorithm>

graph_cluster::graph_cluster(splice_graph &g, vector<hit> &h, int max_gap, bool b)
	: gr(g), hits(h), max_partition_gap(max_gap), store_hits(b)
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
	paired.clear();

	vector<PI> fs;
	build_paired_reads(hits, fs);
	paired.resize(hits.size(), false);

	//printf("total %lu paired reads\n", fs.size());

	for(int i = 0; i < fs.size(); i++)
	{
		int h1 = fs[i].first;
		int h2 = fs[i].second;

		if(hits[h1].rpos > hits[h2].rpos) continue;
		if(hits[h1].pos > hits[h2].pos) continue;

		bool b = consistent_intron_chains(hits[h1].spos, hits[h2].spos);
		if(b == false) continue;

		vector<int> v1;
		vector<int> v2;
		bool b1 = align_hit_to_splice_graph(hits[h1], gr, v1);
		bool b2 = align_hit_to_splice_graph(hits[h2], gr, v2);
		if(b1 == false || b2 == false)  continue;
		if(v1.size() == 0 || v2.size() == 0) continue;

		paired[h1] = true;
		paired[h2] = true;

		PVV pvv(v1, v2);
		if(findex.find(pvv) == findex.end())
		{
			vector<PI> v;
			v.push_back(fs[i]);
			findex.insert(pair<PVV, int>(pvv, groups.size()));
			int32_t p1 = gr.get_vertex_info(v1.front()).lpos;
			int32_t p2 = gr.get_vertex_info(v1.back()).rpos;
			int32_t p3 = gr.get_vertex_info(v2.front()).lpos;
			int32_t p4 = gr.get_vertex_info(v2.back()).rpos;
			groups.push_back(v);
			extend.push_back(p1);
			extend.push_back(p2);
			extend.push_back(p3);
			extend.push_back(p4);
		}
		else
		{
			int k = findex[pvv];
			groups[k].push_back(fs[i]);
		}
	}

	return 0;
}

int graph_cluster::build_pereads_clusters(int g, vector<pereads_cluster> &vc)
{
	const vector<PI> &fs = groups[g];
	vector< vector<int32_t> > vv;
	for(int i = 0; i < fs.size(); i++)
	{
		int h1 = fs[i].first;
		int h2 = fs[i].second;

		vector<int32_t> v;
		v.push_back(hits[h1].pos);
		v.push_back(hits[h1].rpos);
		v.push_back(hits[h2].pos);
		v.push_back(hits[h2].rpos);
		v.push_back(i);
		vv.push_back(v);
	}

	vector< vector<int> > zz = partition(vv, 0);

	for(int i = 0; i < zz.size(); i++)
	{
		if(zz[i].size() == 0) continue;

		int h1 = fs[zz[i][0]].first;
		int h2 = fs[zz[i][0]].second;
		assert(hits[h1].rpos <= hits[h2].rpos);
		assert(hits[h1].pos <= hits[h2].pos);

		bool b = consistent_intron_chains(hits[h1].spos, hits[h2].spos);
		assert(b == true);

		pereads_cluster pc;
		pc.count = 0;
		pc.chain1 = hits[h1].spos;
		pc.chain2 = hits[h2].spos;

		vector<int32_t> bounds(4, 0);
		bounds[0] = hits[h1].pos;
		bounds[1] = hits[h1].rpos;
		bounds[2] = hits[h2].pos;
		bounds[3] = hits[h2].rpos;

		for(int k = 0; k < zz[i].size(); k++)
		{
			h1 = fs[zz[i][k]].first;
			h2 = fs[zz[i][k]].second;

			pc.bounds[0] += hits[h1].pos  - bounds[0];
			pc.bounds[1] += hits[h1].rpos - bounds[1];
			pc.bounds[2] += hits[h2].pos  - bounds[2];
			pc.bounds[3] += hits[h2].rpos - bounds[3];

			pc.count++;
			
			if(store_hits == true)
			{
				pc.hits1.push_back(hits[h1]);
				pc.hits2.push_back(hits[h2]);
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

		/*
		if(pc.bounds[0] == 31579218)
		{
			pc.print(999);
			for(int k = 0; k < zz[i].size(); k++)
			{
				printf("element %d: ", zz[i][k]);
				printv(vv[zz[i][k]]);
				printf("\n");

				h1 = fs[zz[i][k]].first;
				h2 = fs[zz[i][k]].second;
				if(hits[h1].rpos > hits[h2].rpos) continue;
				if(hits[h1].pos > hits[h2].pos) continue;
				hits[h1].print();
				hits[h2].print();
			}
			printf("\n");
		}
		*/
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

vector<bool> graph_cluster::get_paired()
{
	return paired;
}

bool compare_rank0(const vector<int32_t> &x, const vector<int32_t> &y) { return x[0] < y[0]; }
bool compare_rank1(const vector<int32_t> &x, const vector<int32_t> &y) { return x[1] < y[1]; }
bool compare_rank2(const vector<int32_t> &x, const vector<int32_t> &y) { return x[2] < y[2]; }
bool compare_rank3(const vector<int32_t> &x, const vector<int32_t> &y) { return x[3] < y[3]; }
