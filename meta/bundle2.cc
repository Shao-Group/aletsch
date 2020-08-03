/*
Part of aletsch
(c) 2020 by Mingfu Shao, The Pennsylvania State University
See LICENSE for licensing.
*/

#include "bundle2.h"
#include "config.h"
#include "essential.h"

#include <sstream>
#include <algorithm>

bundle2::bundle2(const parameters &c, bundle &&b, int id)
	: cfg(c), bundle(b)
{
	sid = id;
	num_combined = 0;
}

bundle2::bundle2(const parameters &c)
	: cfg(c)
{
	sid = -1;
	num_combined = 0;
}

int bundle2::bridge(const sample_profile &sp)
{
	while(true)
	{
		splice_graph gr;
		graph_builder gb(bb, cfg);
		gb.build(gr);
		gr.build_vertex_index();

		vector<pereads_cluster> vc;
		graph_cluster gc(gr, bb, cfg.max_reads_partition_gap, false);
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
