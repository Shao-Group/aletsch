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

assembler::assembler(const parameters &c, transcript_set_pool &v, transcript_set &tm, mutex &m, thread_pool &p, int r, int g, int i)
	: cfg(c), tspool(v), tmerge(tm), mylock(m), pool(p), rid(r), gid(g), instance(i)
{
	assert(tmerge.rid == rid);
}

int assembler::resolve(vector<bundle*> gv)
{
	int subindex = 0;

	if(gv.size() == 1)
	{
		bundle &bd = *(gv[0]);
		assemble(bd);
	}

	if(gv.size() >= 2)
	{
		// print
		/*
		printf("\n");
		for(int k = 0; k < gv.size(); k++) gv[k]->print(k);

		vector<vector<PID>> sim;
		build_similarity(gv, sim);
		bridge_pairwise(gv, sim);

		sim.clear();
		build_similarity(gv, sim);
		bridge_pairwise(gv, sim);
		*/

		//refine_pairwise(gv, sim);
		bridge(gv);
		//refine(gv);
		assemble(gv);
		//pairwise_assemble(gv, ts, sim, instance);
	}
	return 0;
}

int assembler::build_similarity(vector<bundle*> &gv, vector<vector<PID>> &sim)
{
	vector<vector<int32_t>> splices;
	for(int i = 0; i < gv.size(); i++)
	{
		vector<int32_t> v = gv[i]->hcst.get_splices();
		sort(v.begin(), v.end());
		splices.push_back(std::move(v));
	}

	sim.resize(gv.size());
	for(int i = 0; i < gv.size(); i++)
	{
		for(int j = 0; j < gv.size(); j++)
		{
			if(i == j) continue;

			vector<int32_t> vv(splices[i].size() + splices[j].size(), 0);
			vector<int32_t>::iterator it = set_intersection(splices[i].begin(), splices[i].end(), splices[j].begin(), splices[j].end(), vv.begin());
			int c = it - vv.begin();
			double r = c * 1.0 / splices[i].size();

			if(c <= 0) continue;
			sim[i].push_back(PID(j, c));

			sort(sim[i].begin(), sim[i].end(),
					[](const PID &a, const PID &b) { return a.second > b.second; });
		}
	}
	return 0;
}

int assembler::assemble(bundle &bd)
{
	bd.set_gid(rid, gid, instance, 0);
	splice_graph gr;
	transform(bd, gr, true);

    edge_iterator it;
    PEEI pei = gr.edges();
    for(it = pei.first; it != pei.second; it++)
    {
		edge_descriptor e = (*it);
		int s = e->source();
		int t = e->target();

		//if(s == 0) continue;
		//if(t == gr.num_vertices() - 1) continue;
		//if(gr.get_vertex_info(s).rpos == gr.get_vertex_info(t).lpos) continue;//ignore non-splicing junctions
        
        edge_info ei = gr.get_edge_info(e);
        assert(ei.count == 0);
        assert(ei.samples.size() == 0);
        assert(ei.spAbd.size() == 0);
        ei.samples.insert(bd.sp.sample_id);
        ei.spAbd.insert(make_pair(bd.sp.sample_id, gr.get_edge_weight(e)));
        ei.count = 1;
        gr.set_edge_info(e, ei);
    }
    if(cfg.verbose >= 2)
    {
        printf("\nprint individual graph %s, sample_id=%d\n", gr.gid.c_str(), bd.sp.sample_id);
        gr.print_junction_supports();
    }

	phase_set ps;
	bd.build_phase_set(ps, gr);
	assemble(gr, ps, bd.sp.sample_id);
	bd.clear();
	return 0;
}

int assembler::assemble(vector<bundle*> gv)
{
	assert(gv.size() >= 2);
	int subindex = 0;

	// combined bundle
	bundle bx(cfg, gv[0]->sp);
	bx.copy_meta_information(*(gv[0]));
	for(int k = 0; k < gv.size(); k++) bx.combine(*(gv[k]));
	bx.set_gid(rid, gid, instance, subindex++);

	// combined graph
	splice_graph gx;
	transform(bx, gx, false);	// TODO

	// combined phase set 
	phase_set px;

    map< pair<int32_t, int32_t>, set<int> > junc2sup;
    map< pair< pair<int32_t, int32_t>, int>, double> sup2abd;
    for(int k = 0; k < gv.size(); k++)
    {
        bundle &bd = *(gv[k]);
		splice_graph gr;
		transform(bd, gr, true);

        edge_iterator it;
        PEEI pei = gr.edges();
        for(it = pei.first; it != pei.second; it++)
        {
		    edge_descriptor e = (*it);
		    int s = e->source();
		    int t = e->target();

		    if(s == 0) continue;
		    if(t == gr.num_vertices() - 1) continue;
		    //if(gr.get_vertex_info(s).rpos == gr.get_vertex_info(t).lpos) continue;//ignore non-splicing junctions

            pair<int32_t, int32_t>p = make_pair(gr.get_vertex_info(s).rpos,gr.get_vertex_info(t).lpos);
            junc2sup[p].insert(bd.sp.sample_id);
            
            pair< pair<int32_t, int32_t>, int> psp = make_pair(p, bd.sp.sample_id);
            /*if(sup2abd.find(psp) != sup2abd.end())
                sup2abd[psp] = gr.get_edge_weight(e);
            else*/
            sup2abd[psp] += gr.get_edge_weight(e);
        }
    }

	// assemble individual bundle
	for(int k = 0; k < gv.size(); k++)
	{
		bundle &bd = *(gv[k]);
		bd.set_gid(rid, gid, instance, subindex++);

		splice_graph gr;
		transform(bd, gr, true);

        fix_missing_edges(gr, gx);

        //calculate junction supports based on other samples
        junction_support(gr, junc2sup, sup2abd);
        if(cfg.verbose >= 2) 
        {
            printf("print %d/%lu individual graph %s, sample_id=%d\n", k+1, gv.size(), gr.gid.c_str(), bd.sp.sample_id);
            gr.print_junction_supports();
        }

        for(int j = 0; j < gv.size(); j++)
        {
            bundle &bd1 = *(gv[j]);
            splice_graph gr1;
            transform(bd1, gr1, true);
            start_end_support(bd1.sp.sample_id, gr1, gr);
        }

		phase_set ps;
		bd.build_phase_set(ps, gr);
		px.combine(ps);

		assemble(gr, ps, bd.sp.sample_id);
		//bd.clear();
	}

    for(int k = 0; k < gv.size(); k++)
        gv[k]->clear();

	bx.clear();
	// assemble combined instance
	//assemble(gx, px, ts, -1);
	return 0;
}

int assembler::junction_support(splice_graph &gr, map< pair<int32_t, int32_t>, set<int> > &junc2sup, map< pair< pair<int32_t, int32_t>, int>, double> &sup2abd)
{
    edge_iterator it;
    PEEI pei = gr.edges();
    for(it = pei.first; it != pei.second; it++)
    {
		edge_descriptor e = (*it);
		int s = e->source();
		int t = e->target();

		if(s == 0) continue;
		if(t == gr.num_vertices() - 1) continue;
		//if(gr.get_vertex_info(s).rpos == gr.get_vertex_info(t).lpos) continue;//ignore non-splicing junctions

        pair<int32_t, int32_t>p = make_pair(gr.get_vertex_info(s).rpos,gr.get_vertex_info(t).lpos);
        if(junc2sup.find(p) != junc2sup.end())
        {
            edge_info ei = gr.get_edge_info(e);
            ei.samples = junc2sup[p];
            ei.count = ei.samples.size();
            for(auto sp : ei.samples)
            {
                pair< pair<int32_t, int32_t>, int> psp = make_pair(p, sp);
                if(sup2abd.find(psp) != sup2abd.end()) ei.spAbd[sp] = sup2abd[psp];
            }
            gr.set_edge_info(e, ei);
        }
    }
    return 0;
}

int assembler::start_end_support(int sample_id, splice_graph &gr, splice_graph &gx)
{
    // calculate the support of starting vertices
    edge_iterator it;
	PEEI pei = gr.out_edges(0);
	for(it = pei.first; it != pei.second; it++)
	{
		edge_descriptor e = (*it);
		int s = e->source();
		int t = e->target();
		assert(s == 0);

        //if(gr.edge(t, gr.num_vertices()-1).second) continue;

		int32_t p = gr.get_vertex_info(t).rpos;
		int k = gx.locate_vertex(p - 1);
		if(k < 0) continue;
        int kori = k;

		PEB peb = gx.edge(0, k);
        
        bool cont = true;
        PEB peb2;
        while(!peb.second)
        {
            k--;
            if(k == 0)
            {
                cont = false;
                break;
            }
            if(p - gx.get_vertex_info(k).rpos > 200) cont = false;

            if(gx.get_vertex_info(k+1).lpos != gx.get_vertex_info(k).rpos) cont = false;
            peb2 = gx.edge(k, k+1);
            if(!peb2.second) cont = false;
            if(!cont) break;
            peb = gx.edge(0, k);
        }
        if(!cont) continue;

		edge_info ei = gx.get_edge_info(peb.first);
        ei.samples.insert(sample_id);
        ei.count = ei.samples.size();
        ei.spAbd[sample_id] += gr.get_edge_weight(e);
		gx.set_edge_info(peb.first, ei);
        if(cfg.verbose >= 2) printf("Sample %d supports (%d, %d <- %d(%d))\n", sample_id, 0, k, kori, t);
	}

    // calculate the support of ending vertices
	pei = gr.in_edges(gr.num_vertices() - 1);
	for(it = pei.first; it != pei.second; it++)
	{
		edge_descriptor e = (*it);
		int s = e->source();
		int t = e->target();
		assert(t == gr.num_vertices() - 1);

        //if(gr.edge(0, s).second) continue;

		int32_t p = gr.get_vertex_info(s).lpos;
		//int k = gx.locate_vertex(p - 1)
        int k = gx.locate_vertex(p);
		if(k < 0) continue;
        int kori = k;

		PEB peb = gx.edge(k, gx.num_vertices() - 1);

        bool cont = true;
        PEB peb2;
        while(!peb.second)
        {
            k++;
            if(k == gx.num_vertices() - 1)
            {
                cont = false;
                break;
            }

            if(gx.get_vertex_info(k).lpos - p > 200) cont = false;

            if(gx.get_vertex_info(k-1).rpos != gx.get_vertex_info(k).lpos) cont = false;
            peb2 = gx.edge(k-1, k);
            if(!peb2.second) cont = false;
            if(!cont) break;
            peb = gx.edge(k, gx.num_vertices()-1);
        }
        if(!cont) continue;

		edge_info ei = gx.get_edge_info(peb.first);
        ei.samples.insert(sample_id);
        ei.count = ei.samples.size();
        ei.spAbd[sample_id] += gr.get_edge_weight(e);
		gx.set_edge_info(peb.first, ei);
        if(cfg.verbose >= 2) printf("Sample %d supports (%d(%d) -> %d, %ld)\n", sample_id, kori, s, k, gx.num_vertices()-1);
	}

    return 0;
}

/*int assembler::junction_support(int sample_id, splice_graph &gr, splice_graph &gx)
{
	edge_iterator it;
	PEEI pei = gx.edges();
	for(it = pei.first; it != pei.second; it++)
	{
		edge_descriptor e = (*it);
		int s = e->source();
		int t = e->target();

		if(s == 0) continue;
		if(t == gx.num_vertices() - 1) continue;
		if(gx.get_vertex_info(s).rpos == gx.get_vertex_info(t).lpos) continue;//ignore non-splicing junctions

		int kl = gr.locate_rbound(gx.get_vertex_info(s).rpos);
		int kr = gr.locate_lbound(gx.get_vertex_info(t).lpos);
		if(kl == -1 || kr == -1) continue;

		if(gr.edge(kl, kr).second == false) continue;

		edge_info ei = gx.get_edge_info(e);
        ei.samples.insert(sample_id);
        ei.count = ei.samples.size();
        gx.set_edge_info(e, ei);
    }
    return 0;
}*/

int assembler::transform(bundle &cb, splice_graph &gr, bool revising)
{
	graph_builder gb(cb, cfg, cb.sp);
	gb.build(gr);
	gr.gid = cb.gid;
	gr.build_vertex_index();

	if(revising == true)
	{
		identify_boundaries(gr, cfg);
		remove_false_boundaries(gr, cb, cfg);
		refine_splice_graph(gr);
	}
	return 0;
}

int assembler::fix_missing_edges(splice_graph &gr, splice_graph &gx)
{
	// checking out-edges of 0
	PEEI pi = gr.out_edges(0);
	for(edge_iterator it = pi.first; it != pi.second; it++)
	{
		edge_descriptor e = (*it);
		int t = e->target();
		vertex_info vt = gr.get_vertex_info(t);
		double wt = gr.get_vertex_weight(t);
		int v = gx.locate_rbound(vt.rpos);
		if(v == -1) continue;
		if(gx.in_degree(v) != 1) continue;
		vertex_info vv = gx.get_vertex_info(v);
		edge_descriptor uv = *(gx.in_edges(v).first);
		int u = uv->source();
		double wuv = gx.get_edge_weight(uv);
		if(u == 0) continue;
		vertex_info vu = gx.get_vertex_info(u);
		if(vu.rpos == vv.lpos) continue;
		int s = gr.locate_rbound(vu.rpos);
		if(s == -1) continue;

		int32_t gap = vt.lpos - vv.lpos;
		if(cfg.verbose >= 2) printf("fixing starting boundary t = %d-%d using u = %d-%d, v = %d-%d, gap = %d, wt = %.1lf, wuv = %.1lf\n", 
				vt.lpos, vt.rpos, vu.lpos, vu.rpos, vv.lpos, vv.rpos, gap, wt, wuv);
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
		int unbridged = bd.count_unbridged();

		assert(vc.size() == bs.opt.size());
		for(int j = 0; j < vc.size(); j++)
		{
			if(bs.opt[j].type <= 0) continue;
			cnt1 += 1;
			cnt2 += bd.update_bridges(vc[j].frlist, bs.opt[j].chain, bs.opt[j].strand);
			//vc[j].print(j);
		}

		//if(cfg.verbose >= 2) 
		//printf("gid %s: further bridge %d / %lu clusters, %d / %d fragments\n", bd.gid.c_str(), cnt1, vc.size(), cnt2, unbridged);
	}
	return 0;
}

int assembler::assemble(splice_graph &gx, phase_set &px, int sid)
{
	gx.extend_strands();

	map<int32_t, int32_t> smap, tmap;
	group_start_boundaries(gx, smap, cfg.max_group_boundary_distance);
	group_end_boundaries(gx, tmap, cfg.max_group_boundary_distance);
	px.project_boundaries(smap, tmap);

	hyper_set hx(gx, px);
	hx.filter_nodes(gx);

	if(cfg.verbose >= 2) gx.print();
	if(cfg.verbose >= 2) hx.print_nodes();

	/*
	if(gx.num_vertices() <= 40) 
	{
		string texfile = "tex/" + gx.gid + ".tex";
		gx.draw(texfile);
	}
	*/

	mutex &mt = mylock;	
	transcript_set_pool &tsp = tspool;
	transcript_set &tm = tmerge;
	parameters pa = cfg;
	for(int k = 0; k < cfg.assembly_repeats; k++)
	{
		boost::asio::post(pool, [gx, hx, k, sid, pa, &mt, &tsp, &tm] {

			splice_graph gr(gx);
			hyper_set hs(hx);

			transcript_set ts(gr.chrm, tm.rid, pa.min_single_exon_clustering_overlap);

			//printf("A: tm.rid = %d, ts.rid = %d, this->rid = %d\n", tm.rid, ts.rid, this->rid);

			gr.gid = gx.gid + "." + tostring(k);
			scallop sx(gr, hs, pa, k == 0 ? false : true);
			sx.assemble();

			int z = 0;
			for(int i = 0; i < sx.trsts.size(); i++)
			{
				transcript &t = sx.trsts[i];
				z++;
				t.RPKM = 0;
				ts.add(t, 1, sid, TRANSCRIPT_COUNT_ADD_COVERAGE_ADD);
			}

			if(pa.verbose >= 2) printf("assemble %s: %d transcripts, graph with %lu vertices and %lu edges\n", gr.gid.c_str(), z, gr.num_vertices(), gr.num_edges());
			if(gr.num_vertices() >= 1000) printf("assemble %s: %d transcripts, large graph with %lu vertices and %lu edges\n", gr.gid.c_str(), z, gr.num_vertices(), gr.num_edges());

			mt.lock();
			//tsp.tsets.push_back(ts);
			//printf("B: tm.rid = %d, ts.rid = %d, this->rid = %d\n", tm.rid, ts.rid, this->rid);
			tm.add(ts, TRANSCRIPT_COUNT_ADD_COVERAGE_ADD);
			mt.unlock();
		});
	}

	return 0;
}
