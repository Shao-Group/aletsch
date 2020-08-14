/*
Part of aletsch
(c) 2020 by Mingfu Shao, The Pennsylvania State University
See LICENSE for licensing.
*/

#include "assembler.h"
#include "scallop.h"
#include "graph_builder.h"
#include "graph_cluster.h"
#include "graph_reviser.h"
#include "bridge_solver.h"
#include "essential.h"
#include "constants.h"
#include "filter.h"

#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <thread>
#include <ctime>
#include <boost/asio/post.hpp>
#include <boost/asio/thread_pool.hpp>
#include <boost/pending/disjoint_sets.hpp>

assembler::assembler(const parameters &p)
	: cfg(p)
{
}

int assembler::resolve(vector<bundle*> gv, transcript_set &ts, int instance)
{
	int subindex = 0;

	if(gv.size() == 1)
	{
		bundle &bd = *(gv[0]);
		assemble(bd, ts, instance);
	}

	if(gv.size() >= 2)
	{
		bridge(gv);
		assemble(gv, ts, instance);
	}
	return 0;
}

int assembler::assemble(bundle &bd, transcript_set &ts, int instance)
{
	bd.set_gid(instance, 0);
	splice_graph gr;
	transform(bd, gr, true);

	phase_set ps;
	bd.build_phase_set(ps, gr);
	assemble(gr, ps, ts, bd.sp.sample_id);
	return 0;
}

int assembler::assemble(vector<bundle*> gv, transcript_set &ts, int instance)
{
	assert(gv.size() >= 2);
	int subindex = 0;

	// combined bundle
	bundle bx(cfg, gv[0]->sp);
	bx.copy_meta_information(*(gv[0]));
	for(int k = 0; k < gv.size(); k++) bx.combine(*(gv[k]));
	bx.set_gid(instance, subindex++);

	// combined graph
	splice_graph gx;
	transform(bx, gx, false);	// TODO

	// combined phase set 
	phase_set px;

	// assemble individual bundle
	for(int k = 0; k < gv.size(); k++)
	{
		bundle &bd = *(gv[k]);
		bd.set_gid(instance, subindex++);

		splice_graph gr;
		transform(bd, gr, true);

		phase_set ps;
		bd.build_phase_set(ps, gr);
		px.combine(ps);

		assemble(gr, ps, ts, bd.sp.sample_id);
	}

	// assemble combined instance
	assemble(gx, px, ts, -1);
	return 0;
}

int assembler::transform(bundle &cb, splice_graph &gr, bool revising)
{
	graph_builder gb(cb, cfg, cb.sp);
	if(revising == true) gb.fpe = true;

	gb.build(gr);
	gr.gid = cb.gid;
	gr.build_vertex_index();

	if(revising == true)
	{
		identify_boundaries(gr, cfg);
		remove_false_boundaries(gr, cb);
		refine_splice_graph(gr);
	}
	return 0;
}

int assembler::bridge(vector<bundle*> gv)
{
	assert(gv.size() >= 2);

	// construct combined bundle
	bundle cb(cfg, gv[0]->sp);
	cb.copy_meta_information(*(gv[0]));
	for(int k = 0; k < gv.size(); k++) cb.combine(*(gv[k]));

	// construct combined graph
	splice_graph gr;
	transform(cb, gr, false);

	// bridge each individual bundle
	for(int k = 0; k < gv.size(); k++)
	{
		bundle &bd = *(gv[k]);
		vector<pereads_cluster> vc;
		graph_cluster gc(gr, bd, cfg.max_reads_partition_gap, false);
		gc.build_pereads_clusters(vc);

		if(vc.size() <= 0) continue;

		bridge_solver bs(gr, vc, cfg, bd.sp.insertsize_low, bd.sp.insertsize_high);

		int cnt1 = 0;
		int cnt2 = 0;
		int unbridged = 0;
		for(int j = 0; j < bd.frgs.size(); j++)
		{
			if(bd.frgs[j][2] <= 0) unbridged++;
		}

		assert(vc.size() == bs.opt.size());
		for(int j = 0; j < vc.size(); j++)
		{
			if(bs.opt[j].type <= 0) continue;
			cnt1 += 1;
			cnt2 += bd.update_bridges(vc[j].frlist, bs.opt[j].chain);
		}
		printf("further bridge %d / %lu clusters, %d / %d fragments\n", cnt1, vc.size(), cnt2, unbridged);
	}
	return 0;
}

int assembler::assemble(splice_graph &gx, phase_set &px, transcript_set &ts, int sid)
{
	gx.extend_strands();

	/*
	map<int32_t, int32_t> smap, tmap;
	group_start_boundaries(gx, smap, cfg.max_group_boundary_distance);
	group_end_boundaries(gx, tmap, cfg.max_group_boundary_distance);
	px.project_boundaries(smap, tmap);
	*/

	hyper_set hx(gx, px);
	hx.filter_nodes(gx);

	if(cfg.verbose >= 2) gx.print();

	/*
	if(gx.num_vertices() <= 40) 
	{
		string texfile = "tex/" + gx.gid + ".tex";
		gx.draw(texfile);
	}
	*/

	scallop sx(gx, hx, cfg);
	sx.assemble();

	int z = 0;
	for(int k = 0; k < sx.trsts.size(); k++)
	{
		transcript &t = sx.trsts[k];
		z++;
		t.RPKM = 0;
		ts.add(t, 1, sid, TRANSCRIPT_COUNT_ADD_COVERAGE_ADD);
	}

	if(cfg.verbose >= 2) printf("assemble %s: %d transcripts, graph with %lu vertices and %lu edges, phases = %lu\n", gx.gid.c_str(), z, gx.num_vertices(), gx.num_edges(), px.pmap.size());
	//gx.print();

	return 0;
}
