#include "graph_hits.h"
#include "essential.h"
#include "util.h"

#include <algorithm>

graph_hits::graph_hits(splice_graph &g, vector<hit> &h)
	: gr(g), hits(h)
{} 

int graph_hits::build_paired_reads_clusters(vector<PRC> &vpr, vector<bool> &paired)
{
	typedef pair< vector<int>, vector<int> > PVV;
	map<PVV, int> findex;			// index for fclusters

	vpr.clear();
	paired.assign(hits.size(), false);

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

		int32_t k1l = hits[h1].pos - gr.get_vertex_info(v1.front()).lpos;
		int32_t k1r = gr.get_vertex_info(v1.back()).rpos - hits[h1].rpos;
		int32_t k2l = hits[h2].pos - gr.get_vertex_info(v2.front()).lpos;
		int32_t k2r = gr.get_vertex_info(v2.back()).rpos - hits[h2].rpos;

		PVV pvv(v1, v2);
		if(findex.find(pvv) == findex.end())
		{
			rcluster r1;
			rcluster r2;
			r1.vv = v1;
			r2.vv = v2;
			r1.vl.push_back(k1l);
			r1.vr.push_back(k1r);
			r2.vl.push_back(k2l);
			r2.vr.push_back(k2r);

			findex.insert(pair<PVV, int>(pvv, vpr.size()));
			vpr.push_back(PRC(r1, r2));
		}
		else
		{
			int k = findex[pvv];
			assert(vpr[k].first.vv == v1);
			assert(vpr[k].second.vv == v2);
			vpr[k].first.vl.push_back(k1l);
			vpr[k].first.vr.push_back(k1r);
			vpr[k].second.vl.push_back(k2l);
			vpr[k].second.vr.push_back(k2r);
		}
	}

	return 0;
}

int graph_hits::build_hyper_set_from_unpaired_reads(const vector<bool> &paired, hyper_set &hs)
{
	for(int i = 0; i < hits.size(); i++)
	{
		if(paired[i] == true) continue;
		vector<int> v;
		bool b = align_hit_to_splice_graph(hits[i], gr, v);
		if(b == false) continue;
		for(int k = 0; k < v.size(); k++) v[k]--;
		hs.add_node_list(v, 1);
	}
	return 0;
}
