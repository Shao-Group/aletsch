/*
Part of aletsch
(c) 2020 by Mingfu Shao, The Pennsylvania State University
See LICENSE for licensing.
*/

#include "bundle.h"
#include "config.h"
#include "essential.h"
#include "graph_builder.h"
#include "graph_cluster.h"
#include "bridge_solver.h"

#include <sstream>
#include <algorithm>

bundle::bundle(const parameters &c, const sample_profile &s)
	: cfg(c), sp(s)
{
	num_combined = 0;
}

bundle::bundle(const parameters &c, const sample_profile &s, bundle_base &&bb)
	: cfg(c), sp(s), bundle_base(bb)
{
	num_combined = 0;
}

int bundle::set_gid(int instance, int subindex)
{
	char name[10240];
	sprintf(name, "instance.%d.%d", instance, subindex);
	gid = name;
	return 0;
}

int bundle::set_gid(int rid, int g, int instance, int subindex)
{
	char name[10240];
	sprintf(name, "instance.%d.%d.%d.%d", rid, g, instance, subindex);
	gid = name;
	return 0;
}

int bundle::copy_meta_information(const bundle &bb)
{
	chrm = bb.chrm;
	strand = bb.strand;
	tid = bb.tid;
	lpos = bb.lpos;
	rpos = bb.rpos;
	return 0;
}

int bundle::bridge()
{
	while(true)
	{
		splice_graph gr;
		graph_builder gb(*this, cfg, sp);
		gb.build(gr);
		gr.build_vertex_index();

		vector<pereads_cluster> vc;
		graph_cluster gc(gr, *this, cfg.max_reads_partition_gap, false);
		gc.build_pereads_clusters(vc);

		bridge_solver bs(gr, vc, cfg, sp.insertsize_low, sp.insertsize_high);

		int cnt = 0;
		assert(vc.size() == bs.opt.size());
		for(int k = 0; k < vc.size(); k++)
		{
			if(bs.opt[k].type <= 0) continue;
			cnt += update_bridges(vc[k].frlist, bs.opt[k].chain);
		}

		//printf("total frags %lu, bridged frags = %d\n", bb.frgs.size(), cnt);
		if(cnt <= 0) break;
	}
	return 0;
}

int bundle::combine(const bundle &bb)
{
	num_combined += bb.num_combined;
	assert(strand == bb.strand);
	assert(chrm == bb.chrm);
	assert(tid == bb.tid);
	if(lpos > bb.lpos) lpos = bb.lpos;
	if(rpos < bb.rpos) rpos = bb.rpos;
	hcst.add(bb.hcst);
	fcst.add(bb.fcst);
	for(SIMI z = bb.mmap.begin(); z != bb.mmap.end(); z++) mmap += *z;
	for(SIMI z = bb.imap.begin(); z != bb.imap.end(); z++) imap += *z;
	return 0;
}

int bundle::print(int index)
{
	printf("bundle %d: sample = %d, ", index, sp.sample_id);

	// statistic xs
	int n0 = 0, np = 0, nq = 0;
	for(int i = 0; i < hits.size(); i++)
	{
		if(hits[i].xs == '.') n0++;
		if(hits[i].xs == '+') np++;
		if(hits[i].xs == '-') nq++;
	}

	printf("tid = %d, range = %s:%d-%d, orient = %c, #hits = %lu, +/-/. = %d / %d / %d\n",
			tid, chrm.c_str(), lpos, rpos, strand, hits.size(), np, nq, n0);

	return 0;
}

int bundle::save_borrowed_path(const vector<int32_t> &v, double w)
{
	borrowed_paths.push_back(v);
	borrowed_weights.push_back(w);
	return 0;
}

int bundle::print_borrowed_paths()
{
	for(int k = 0; k < borrowed_paths.size(); k++)
	{
		printf("borrowed path %d: weight = %.2lf, path = ", k, borrowed_weights[k]);
		printv(borrowed_paths[k]);
		printf("\n");
	}
	return 0;
}

int bundle::digest_borrowed_paths()
{
	map<vector<int32_t>, double> m;
	for(int k = 0; k < borrowed_paths.size(); k++)
	{
		vector<int32_t> &v = borrowed_paths[k];
		double w = borrowed_weights[k];
		if(m.find(v) == m.end()) m.insert(make_pair(v, w));
		else m[v] += w;
	}

	for(auto &x: m)
	{
		add_borrowed_path(x.first, x.second);
		printf("add borrowed path, w = %.2lf, p = ", x.second);
		printv(x.first);
		printf("\n");
	}
	return 0;
}

/*
int bundle::count_unbridged_fragments()
{
	int c = 0;
	for(int j = 0; j < frgs.size(); j++)
	{
		if(frgs[j][2] <= 0) c++;
	}
	return c;
}
*/
