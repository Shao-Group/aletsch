#include "graph_hits.h"
#include "essential.h"
#include "util.h"

#include <algorithm>

graph_hits::graph_hits(splice_graph &g, vector<hit> &h)
	: gr(g), hits(h)
{} 

int graph_hits::group_pereads(vector< vector<PI> > &vc, vector<bool> &paired)
{
	typedef pair< vector<int>, vector<int> > PVV;
	map<PVV, int> findex;

	vc.clear();
	paired.clear();

	vector<PI> fs;
	build_paired_reads(hits, fs);

	for(int i = 0; i < fs.size(); i++)
	{
		int h1 = fs[i].first;
		int h2 = fs[i].second;
		vector<int> v1;
		vector<int> v2;
		bool b1 = align_hit_to_splice_graph(hits[h1], gr, v1);
		bool b2 = align_hit_to_splice_graph(hits[h2], gr, v2);

		if(b1 == false || b2 == false)  continue;

		paired[h1] = true;
		paired[h2] = true;

		PVV pvv(v1, v2);
		if(findex.find(pvv) == findex.end())
		{
			vector<PI> v;
			v.push_back(fs[i]);
			findex.insert(pair<PVV, int>(pvv, vc.size()));
			vc.push_back(v);
		}
		else
		{
			int k = findex[pvv];
			vc[k].push_back(fs[i]);
		}
	}

	return 0;
}

int graph_hits::build_pereads_clusters(const vector<PI> &fs, vector<pereads_cluster> &vc)
{
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

	vc.clear();
	for(int i = 0; i < zz.size(); i++)
	{
		if(zz[i].size() == 0) continue;

		int h1 = fs[zz[i][0]].first;
		int h2 = fs[zz[i][0]].second;

		pereads_cluster pc;
		pc.chain1 = hits[h1].spos;
		pc.chain2 = hits[h2].spos;
		pc.count = zz[i].size();

		for(int k = 1; k < zz[i].size(); k++)
		{
			h1 = fs[zz[i][k]].first;
			h2 = fs[zz[i][k]].second;
			pc.bounds[0] += hits[h1].pos;
			pc.bounds[1] += hits[h1].rpos;
			pc.bounds[2] += hits[h2].pos;
			pc.bounds[3] += hits[h2].rpos;
		}

		pc.bounds[0] /= pc.count; 
		pc.bounds[1] /= pc.count;
		pc.bounds[2] /= pc.count; 
		pc.bounds[3] /= pc.count;
	}

	return 0;
}

vector< vector<int> > graph_hits::partition(vector< vector<int32_t> > &fs, int r)
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

	int32_t max_partition_gap = 20;

	int pre = 0;
	for(int k = 1; k <= fs.size(); k++)
	{
		if(k < fs.size() && fs[k][r] - fs[k - 1][r] <= max_partition_gap) continue;

		vector< vector<int32_t> > fs1;
		for(int i = pre; i < k; i++) fs1.push_back(fs[i]);

		vector< vector<int> > vv1 = partition(fs1, r + 1);
		vv.insert(vv.end(), vv1.begin(), vv1.end());

		pre = k;
	}
	return vv;
}

int graph_hits::build_phase_set_from_unpaired_reads(const vector<bool> &paired, phase_set &ps)
{
	for(int i = 0; i < hits.size(); i++)
	{
		if(paired[i] == true) continue;
		vector<int> v;
		bool b = align_hit_to_splice_graph(hits[i], gr, v);
		if(b == false) continue;
		vector<int32_t> p;
		build_exon_coordinates_from_path(gr, v, p);
		ps.add(p, 1);
	}
	return 0;
}

bool compare_rank0(const vector<int32_t> &x, const vector<int32_t> &y) { return x[0] < y[0]; }
bool compare_rank1(const vector<int32_t> &x, const vector<int32_t> &y) { return x[1] < y[1]; }
bool compare_rank2(const vector<int32_t> &x, const vector<int32_t> &y) { return x[2] < y[2]; }
bool compare_rank3(const vector<int32_t> &x, const vector<int32_t> &y) { return x[3] < y[3]; }
