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

int assembler::assemble(vector<bundle*> gv, int batch, int instance, transcript_set &ts)
{
	int subindex = 0;

	if(gv.size() == 1)
	{
		bundle &gt = *(gv[0]);
		gt.set_gid(batch, instance, subindex++);

		// TODO
		//if(cfg.boost_precision == true) gt.refine_junctions();

		assemble(gt, ts, TRANSCRIPT_COUNT_ADD_COVERAGE_ADD);
		//ts.increase_count(1);
		// TODO, write reads
	}
	else if(gv.size() >= 2)
	{
		bundle cx(cfg, gv[0]->sp);
		resolve_cluster(gv, cx);
		cx.set_gid(batch, instance, subindex++);
		assemble(cx, ts, TRANSCRIPT_COUNT_ADD_COVERAGE_ADD);

		for(int i = 0; i < gv.size(); i++)
		{
			gv[i]->set_gid(batch, instance, subindex++);
			assemble(*gv[i], ts, TRANSCRIPT_COUNT_ADD_COVERAGE_ADD);
		}
	}

	return 0;
}

int assembler::assemble(bundle &cb, transcript_set &ts, int mode)
{
	splice_graph gr;
	graph_builder gb(cb, cfg);
	gb.build(gr);
	gr.gid = cb.gid;

	vector<transcript> vt;

	// TODO build phase-set
	//assemble(gr, cb.ps, vt, cb.num_combined);

	for(int k = 0; k < vt.size(); k++)
	{
		ts.add(vt[k], 1, cb.sp.sample_id, mode);
	}
	return 0;
}

int assembler::assemble(splice_graph &gx, phase_set &px, vector<transcript> &vt, int combined)
{
	gx.build_vertex_index();
	gx.extend_strands();

	map<int32_t, int32_t> smap, tmap;
	group_start_boundaries(gx, smap, cfg.max_group_boundary_distance);
	group_end_boundaries(gx, tmap, cfg.max_group_boundary_distance);
	px.project_boundaries(smap, tmap);

	refine_splice_graph(gx);
	hyper_set hx(gx, px);
	hx.filter_nodes(gx);

	/*
	gx.print();
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
		vt.push_back(t);
	}

	printf("assemble %s: %d transcripts, combined = %d, graph with %lu vertices and %lu edges, phases = %lu\n", gx.gid.c_str(), z, combined, gx.num_vertices(), gx.num_edges(), px.pmap.size());
	gx.print();

	/*
	for(int k = 0; k < sx.trsts.size(); k++)
	{
		transcript &t = sx.trsts[k];
		t.write(cout);
	}
	printf("\n");
	*/

	return 0;
}

int assembler::resolve_cluster(vector<bundle*> gv, bundle &cb)
{
	assert(gv.size() >= 2);

	// construct combined bundle
	cb.copy_meta_information(*(gv[0]));
	for(int k = 0; k < gv.size(); k++) cb.combine(*(gv[k]));

	// construct combined graph
	splice_graph gr;
	graph_builder gb(cb, cfg);
	gb.build(gr);
	gr.build_vertex_index();

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

	// TODO write reads
	return 0;
}
