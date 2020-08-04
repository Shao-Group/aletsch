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

bundle::bundle(const parameters &c, sample_profile &s, bundle_base &&bb)
	: cfg(c), sp(s), bundle_base(bb)
{
	num_combined = 0;
}

int bundle::set_gid(int batch, int instance, int subindex)
{
	char name[10240];
	sprintf(name, "instance.%d.%d.%d", batch, instance, subindex);
	gid = name;
	return 0;
}

int bundle::bridge()
{
	while(true)
	{
		splice_graph gr;
		graph_builder gb(*this, cfg);
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
	assert(chrm == bb.chrm);
	return 0;
}
