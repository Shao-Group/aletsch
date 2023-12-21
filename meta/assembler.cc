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

	// FEATURE: gv.size()
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

    gr.reads = bd.frgs.size();
    gr.subgraph = 1;


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
        
        edge_info &ei = gr.get_editable_edge_info(e);
        assert(ei.count == 0);
        assert(ei.samples.size() == 0);
        assert(ei.spAbd.size() == 0);
        ei.samples.insert(bd.sp.sample_id);
        ei.spAbd.insert(make_pair(bd.sp.sample_id, gr.get_edge_weight(e)));
        ei.abd = gr.get_edge_weight(e);
        ei.count = 1;
        //gr.set_edge_info(e, ei);
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

    gx.reads = bx.frgs.size();
    gx.subgraph = gv.size();
	//printf("Merged #reads: hits = %lu, frgs = %lu, gx.reads = %d\n", bx.hits.size(), bx.frgs.size(), gx.reads);

	// combined phase set 
	phase_set px;

    //junction supports and supported sample abundance
    unordered_map<int64_t, set<int> > junc2sup;
    unordered_map<int64_t, unordered_map<int, double>> sup2abd;

    // combined support
    edge_iterator itx;
    PEEI peix = gx.edges();
    for(itx = peix.first; itx != peix.second; itx++)
    {
        edge_descriptor e = (*itx);
        int s = e->source();
        int t = e->target();

        edge_info & ei = gx.get_editable_edge_info(e);
        ei.samples.clear();
        ei.spAbd.clear();
        ei.samples.insert(-1);
        ei.spAbd.insert(make_pair(-1, gx.get_edge_weight(e)));
        //printf("bx.sp.sample_id:%d\n", bx.sp.sample_id);
        ei.abd = gx.get_edge_weight(e);
        ei.count = 1;

        if(s == 0) continue;
        if(t == gx.num_vertices() - 1) continue;

        pair<int32_t, int32_t>p0 = make_pair(gx.get_vertex_info(s).rpos, gx.get_vertex_info(t).lpos);
        if(p0.first == p0.second) continue;//ignore non-splicing junctions

		int64_t p = pack(p0.first, p0.second);
        junc2sup[p].insert(-1);
        
        //pair<int64_t, int> psp = make_pair(p, -1);
        //sup2abd[psp] += gx.get_edge_weight(e);
        sup2abd[p].insert(make_pair(-1, gx.get_edge_weight(e)));
    }

    //transform individual bundle to individual graph
    vector<splice_graph*> grv;
    vector<int> idv;
    size_t max_v_num = 0;

    // individual junction supports
    for(int k = 0; k < gv.size(); k++)
    {
        bundle &bd = *(gv[k]);
        bd.set_gid(rid, gid, instance, subindex++);

		splice_graph* grp = new splice_graph();
        grv.push_back(grp);
        splice_graph& gr = *grp; 
        transform(bd, gr, true);
        idv.push_back(bd.sp.sample_id);

        gr.reads = bd.frgs.size();
        gr.subgraph = gv.size();
        max_v_num = max(max_v_num, gr.num_vertices());
        //printf("Graph %d, #reads: hits = %lu, frgs = %lu, gx.reads = %d\n", k+1, bd.hits.size(), bd.frgs.size(), gr.reads);

        edge_iterator it;
        PEEI pei = gr.edges();
        for(it = pei.first; it != pei.second; it++)
        {
		    edge_descriptor e = (*it);
		    int s = e->source();
		    int t = e->target();

            edge_info & ei = gr.get_editable_edge_info(e);
            ei.samples.clear();
            ei.spAbd.clear();
            ei.samples.insert(bd.sp.sample_id);
            ei.spAbd.insert(make_pair(bd.sp.sample_id, gr.get_edge_weight(e)));
            ei.abd = gr.get_edge_weight(e);
            ei.count = 1;

            if(s == 0) continue;
		    if(t == gr.num_vertices() - 1) continue;

            pair<int32_t, int32_t> p0 = make_pair(gr.get_vertex_info(s).rpos,gr.get_vertex_info(t).lpos);
            if(p0.first == p0.second) continue;//ignore non-splicing junctions
			int64_t p = pack(p0.first, p0.second);
            junc2sup[p].insert(bd.sp.sample_id);
            
            //pair<int64_t, int> psp = make_pair(p, bd.sp.sample_id);
            //sup2abd[psp] += gr.get_edge_weight(e);
            sup2abd[p].insert(make_pair(bd.sp.sample_id, gr.get_edge_weight(e)));
        }
    }

    //append gx to grv
    grv.push_back(&gx);
    idv.push_back(-1);
    
    //start_end_support(grv, idv);

    //assemble merged graph when the largest graph in the bundle has <=150 vertices
    bool assemble_merged = true;
    //if(max_v_num <= 150) assemble_merged = true;

	// assemble individual bundle
	for(int k = 0; k < gv.size(); k++)
	{
        bundle &bd = *(gv[k]);
		//bd.set_gid(rid, gid, instance, subindex++);
        //splice_graph gr;
		//transform(bd, gr, true);
        splice_graph &gr = *(grv[k]);


        fix_missing_edges(gr, gx);

		if(cfg.verbose >= 2) 
        {
            printf("\n-----preprocess splice graph %s, vertices = %lu, edges = %lu\n", gr.gid.c_str(), gr.num_vertices(), gr.num_edges());
            gr.print();
        }
        //calculate junction supports based on other samples
        junction_support(gr, junc2sup, sup2abd);
        for(int j = 0; j < gv.size(); j++)
        {
            bundle &bd1 = *(gv[j]);
            //splice_graph gr1;
            //transform(bd1, gr1, true);
            splice_graph &gr1 = *(grv[j]);

            start_end_support(bd1.sp.sample_id, gr1, gr);
            non_splicing_support(bd1.sp.sample_id, gr1, gr);
            boundary_extend(bd1.sp.sample_id, gr, gr1, 1);
            boundary_extend(bd1.sp.sample_id, gr, gr1, 2);
            boundary_extend(bd1.sp.sample_id, gr, gr1, 3);
        }

        if(cfg.verbose >= 2) 
        {
            printf("print %d/%lu individual graph %s, sample_id=%d\n", k+1, gv.size(), gr.gid.c_str(), bd.sp.sample_id);
            gr.print_junction_supports();
        }
		phase_set ps;
		bd.build_phase_set(ps, gr);
		px.combine(ps);

        //calculate start&end&non-splicing suppots for combined graph
        if(assemble_merged)
        {
            start_end_support(bd.sp.sample_id, gr, gx);
            non_splicing_support(bd.sp.sample_id, gr, gx);
            boundary_extend(-1, gr, gx, 1);
        }

		assemble(gr, ps, bd.sp.sample_id);
		//bd.clear();
	}

    for(int k = 0; k < gv.size(); k++)
    {
        gv[k]->clear();
        delete grv[k];
    }

	bx.clear();

    if(assemble_merged)
    {
        junction_support(gx, junc2sup, sup2abd);
        //start_end_support(-1, gx, gx);
        //non_splicing_support(-1, gx, gx);

        if(cfg.verbose >= 2) 
        {
            printf("print combined graph %s\n", gx.gid.c_str());
            gx.print_junction_supports();
        }

        // assemble combined instance
        assemble(gx, px, -1);
    }
    return 0;
}

int assembler::junction_support(splice_graph &gr, unordered_map<int64_t, set<int> > &junc2sup, unordered_map<int64_t, unordered_map<int, double>> &sup2abd)
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
		if(gr.get_vertex_info(s).rpos == gr.get_vertex_info(t).lpos) continue;//ignore non-splicing junctions

        //pair<int32_t, int32_t>p = make_pair(gr.get_vertex_info(s).rpos,gr.get_vertex_info(t).lpos);
		int64_t p = pack(gr.get_vertex_info(s).rpos,gr.get_vertex_info(t).lpos);
        if(junc2sup.find(p) != junc2sup.end())
        {
            edge_info &ei = gr.get_editable_edge_info(e);
            ei.samples = junc2sup[p];
			ei.spAbd = sup2abd[p];
            ei.count = ei.samples.size();
			for(auto &z : sup2abd[p])
			{
				ei.abd += z.second;
			}
			/*
            for(auto sp : ei.samples)
            {
                pair<int64_t, int> psp = make_pair(p, sp);
                if(sup2abd.find(psp) != sup2abd.end()) 
                {
                    //assert(ei.spAbd[sp] < SMIN);
                    ei.spAbd[sp] = sup2abd[psp];
                    ei.abd += ei.spAbd[sp]; 
                }
            }
			*/
            //gr.set_edge_info(e, ei);
        }
    }
    return 0;
}

int assembler::non_splicing_support(int sample_id, splice_graph &gr, splice_graph &gx)
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
        edge_info &ei = gx.get_editable_edge_info(e);

		if(gx.get_vertex_info(s).rpos == gx.get_vertex_info(t).lpos)
        {
           int32_t p = gx.get_vertex_info(t).lpos;
           int k1 = gr.locate_vertex(p-1);
           int k2 = gr.locate_vertex(p);
           if(k1 < 0 || k2 < 0) continue;

           if(k1 == k2)
           {
               ei.samples.insert(sample_id);
               ei.count = ei.samples.size();
               ei.spAbd[sample_id] += gr.get_vertex_weight(k1);
               ei.abd += gr.get_vertex_weight(k1);
               //gx.set_edge_info(e, ei);
               if(cfg.verbose >= 3) printf("Non-splicing edge(%d, %d) supported by vertex %d, sample_id=%d, weight=%.2f\n", s, t, k1, sample_id, gr.get_vertex_weight(k1));
           }
           else if(gr.get_vertex_info(k1).rpos==gr.get_vertex_info(k2).lpos && gr.edge(k1, k2).second)
           {
               ei.samples.insert(sample_id);
               ei.count = ei.samples.size();
               ei.spAbd[sample_id] += gr.get_edge_weight(gr.edge(k1, k2).first);
               ei.abd += gr.get_edge_weight(gr.edge(k1, k2).first);
               //gx.set_edge_info(e, ei);
               if(cfg.verbose >= 3) printf("Non-splicing edge(%d, %d) supported by edge(%d, %d), sample_id=%d, weight=%.2f\n", s, t, k1, k2, sample_id, gr.get_edge_weight(gr.edge(k1, k2).first));
           }

        }
    }
    return 0;
}

int assembler::start_end_support(vector<splice_graph*> &grv, const vector<int> &idv)
{
	// build start interval_set_map
	interval_set_pair_map ism_start;
	int32_t gap = 200;
	for(int k = 0; k < grv.size(); k++)
	{
		splice_graph &gr = *(grv[k]);
		set<PI> se;

		edge_iterator it;
		PEEI pei = gr.out_edges(0);
		for(it = pei.first; it != pei.second; it++)
		{
			edge_descriptor e = (*it);
			int s = e->source();
			int t = e->target();
			assert(s == 0);
			se.insert(PI(k, t));

			int32_t p1 = gr.get_vertex_info(t).lpos;
			int32_t p2 = gr.get_vertex_info(t).rpos;

			while(p2 - p1 < gap)
			{
				if(t + 1 >= gr.num_vertices() - 1) break;	
				PEB peb = gr.edge(t, t + 1);
				if(peb.second == false) break;
				if(gr.get_vertex_info(t + 1).lpos != gr.get_vertex_info(t).rpos) break;
				p2 = gr.get_vertex_info(t + 1).rpos;
				t++;
			}

			if(p2 > p1 + gap) p2 = p1 + gap;
			ism_start += make_pair(interval32(p1, p2), se);
		}
	}

	// first collect all pairs
	map<PI, vector<PI>> mpairs;

	for(ISPMI it = ism_start.begin(); it != ism_start.end(); it++)
	{
		vector<PI> v1(it->second.size());
		vector<PI> v2(it->second.size());

		if(it != ism_start.begin())
		{
			auto ix = it;
			--ix;

			auto i1 = std::set_difference(it->second.begin(), it->second.end(), ix->second.begin(), ix->second.end(), v1.begin());
			v1.resize(i1 - v1.begin());

			auto i2 = std::set_difference(it->second.begin(), it->second.end(), v1.begin(), v1.end(), v2.begin());
			v2.resize(i2 - v2.begin());

			assert(v1.size() + v2.size() == it->second.size());

			for(int i = 0; i < v1.size(); i++)
			{
				assert(mpairs.find(v1[i]) == mpairs.end());
				vector<PI> v;
				v.insert(v.end(), v1.begin() + i + 1, v1.end());
				v.insert(v.end(), v2.begin(), v2.end());
				mpairs.insert(make_pair(v1[i], v));
			}
		}
	}

	for(auto x: mpairs)
	{
		auto &p = x.first;
		auto &v = x.second;

		int si = idv[p.first];
		int starti = p.second;
		splice_graph &gi = *(grv[p.first]);
		PEB pi = gi.edge(0, starti);
		assert(pi.second == true);
		edge_info &ei = gi.get_editable_edge_info(pi.first);

		for(int j = 0; j < v.size(); j++)
		{
			int sj = idv[v[j].first];
			if(si == sj) continue;

			int startj = v[j].second;
			splice_graph &gj = *(grv[v[j].first]);
			PEB pj = gj.edge(0, startj);
			assert(pj.second == true);
			edge_info &ej = gj.get_editable_edge_info(pj.first);

			ei.samples.insert(sj);
			ej.samples.insert(si);

			ei.count = ei.samples.size();
			ej.count = ej.samples.size();

			ei.spAbd[sj] += gj.get_edge_weight(pj.first);
			ej.spAbd[si] += gi.get_edge_weight(pi.first);

			ei.abd += gj.get_edge_weight(pj.first);
			ej.abd += gi.get_edge_weight(pi.first);
		}
	}

    // build end interval_set_map
	interval_set_pair_map ism_end;
	for(int k = 0; k < grv.size(); k++)
	{
		splice_graph &gr = *(grv[k]);
		set<PI> se;

		edge_iterator it;
		PEEI pei = gr.in_edges(gr.num_vertices()-1);
		for(it = pei.first; it != pei.second; it++)
		{
			edge_descriptor e = (*it);
			int s = e->source();
			int t = e->target();
			assert(t == gr.num_vertices()-1);
			se.insert(PI(k, s));

			int32_t p1 = gr.get_vertex_info(s).lpos;
			int32_t p2 = gr.get_vertex_info(s).rpos;

			while(p2 - p1 < gap)
			{
				if(s - 1 <= 0) break;	
				PEB peb = gr.edge(s - 1, s);
				if(peb.second == false) break;
				if(gr.get_vertex_info(s - 1).rpos != gr.get_vertex_info(s).lpos) break;
				p1 = gr.get_vertex_info(s - 1).lpos;
				s--;
			}

			if(p2 > p1 + gap) p1 = p2 - gap;
			ism_end += make_pair(interval32(p1, p2), se);
		}
	}

	// collect all pairs
	mpairs.clear();

	for(ISPMI it = ism_end.begin(); it != ism_end.end(); it++)
	{
		vector<PI> v1(it->second.size());
		vector<PI> v2(it->second.size());

		if(it != ism_end.begin())
		{
			auto ix = it;
			--ix;

			auto i1 = std::set_difference(it->second.begin(), it->second.end(), ix->second.begin(), ix->second.end(), v1.begin());
			v1.resize(i1 - v1.begin());

			auto i2 = std::set_difference(it->second.begin(), it->second.end(), v1.begin(), v1.end(), v2.begin());
			v2.resize(i2 - v2.begin());

			assert(v1.size() + v2.size() == it->second.size());

			for(int i = 0; i < v1.size(); i++)
			{
				assert(mpairs.find(v1[i]) == mpairs.end());
				vector<PI> v;
				v.insert(v.end(), v1.begin() + i + 1, v1.end());
				v.insert(v.end(), v2.begin(), v2.end());
				mpairs.insert(make_pair(v1[i], v));
			}
		}
	}

	for(auto x: mpairs)
	{
		auto &p = x.first;
		auto &v = x.second;

		int si = idv[p.first];
		int endi = p.second;
		splice_graph &gi = *(grv[p.first]);
		PEB pi = gi.edge(endi, gi.num_vertices()-1);
		assert(pi.second == true);
		edge_info &ei = gi.get_editable_edge_info(pi.first);

		for(int j = 0; j < v.size(); j++)
		{
			int sj = idv[v[j].first];
			if(si == sj) continue;

			int endj = v[j].second;
			splice_graph &gj = *(grv[v[j].first]);
			PEB pj = gj.edge(endj, gj.num_vertices()-1);
			assert(pj.second == true);
			edge_info &ej = gj.get_editable_edge_info(pj.first);

			ei.samples.insert(sj);
			ej.samples.insert(si);

			ei.count = ei.samples.size();
			ej.count = ej.samples.size();

			ei.spAbd[sj] += gj.get_edge_weight(pj.first);
			ej.spAbd[si] += gi.get_edge_weight(pi.first);

			ei.abd += gj.get_edge_weight(pj.first);
			ej.abd += gi.get_edge_weight(pi.first);
		}
	}

    return 0;
}

int assembler::start_end_support(int sample_id, splice_graph &gr, splice_graph &gx)
{
    // calculate support of starting vertices
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

		edge_info &ei = gx.get_editable_edge_info(peb.first);
        ei.samples.insert(sample_id);
        ei.count = ei.samples.size();
        ei.spAbd[sample_id] += gr.get_edge_weight(e);
        ei.abd += gr.get_edge_weight(e);
		//gx.set_edge_info(peb.first, ei);
        if(cfg.verbose >= 3) printf("Sample %d supports (%d, %d <- %d(%d))\n", sample_id, 0, k, kori, t);
	}

    // calculate support of ending vertices
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

		edge_info &ei = gx.get_editable_edge_info(peb.first);
        ei.samples.insert(sample_id);
        ei.count = ei.samples.size();
        ei.spAbd[sample_id] += gr.get_edge_weight(e);
        ei.abd += gr.get_edge_weight(e);
		//gx.set_edge_info(peb.first, ei);
        if(cfg.verbose >= 3) printf("Sample %d supports (%d(%d) -> %d, %ld)\n", sample_id, kori, s, k, gx.num_vertices()-1);
	}

    return 0;
}

//record broken boundaries features of gr, inferred from gx
int assembler::boundary_extend(int sample_id, splice_graph &gr, splice_graph &gx, int pos_type)
{
    //loss of start 
    edge_iterator it;
    PEEI pei = gr.out_edges(0);
    for(it = pei.first; it != pei.second; it++)
    {
        edge_descriptor e = *it;
        int s = e->source();
        int t = e->target();
        assert(s == 0);

        vertex_info &vi = gr.get_editable_vertex_info(t);
        int k = -1;
        if(pos_type == 1)
            k = gx.locate_vertex(vi.lpos);
        else if (pos_type == 2)
            k = gx.locate_vertex(vi.rpos-1);
        else if(pos_type == 3 && gr.edge(t, t+1).second && (gr.get_vertex_info(t).rpos==gr.get_vertex_info(t+1).lpos) && (t+1 < gr.num_vertices()-1))
        {
            k = gx.locate_vertex(vi.rpos);
        }
        if(k <= 0 || gx.edge(0, k).second) continue;
        double new_loss = 0;
        if(gx.edge(k-1, k).second && (gx.get_vertex_info(k-1).rpos == gx.get_vertex_info(k).lpos))
        {
            new_loss = gx.get_in_weights(k)-gx.get_edge_weight(gx.edge(k-1, k).first);
        }
        else
        {
            new_loss = gx.get_in_weights(k);
        }
        if(cfg.verbose >= 2) printf("Start vertex %d(gx=%d, vertex=%d) boundary_loss = %.2lf\n", t, sample_id, k, new_loss);
        if(sample_id == -1 && pos_type == 1)
            vi.boundary_merged_loss += new_loss; 
        else
        {
            if(pos_type == 1)
                vi.boundary_loss1 += new_loss; 
            else if(pos_type == 2)
                vi.boundary_loss2 += new_loss; 
            else if(pos_type == 3)
                vi.boundary_loss3 += new_loss; 
        }
        //gr.set_vertex_info(t, vi);
    }

    //loss of end
    pei = gr.in_edges(gr.num_vertices()-1);
    for(it = pei.first; it != pei.second; it++)
    {
        edge_descriptor e = *it;
        int s = e->source();
        int t = e->target();
        assert(t = gr.num_vertices()-1);

        vertex_info &vi = gr.get_editable_vertex_info(s);
        int k = -1;
        if(pos_type == 1)
            k = gx.locate_vertex(vi.rpos -1);
        else if(pos_type == 2)
            k = gx.locate_vertex(vi.lpos);
        else if(pos_type == 3 && s>1 && gr.edge(s-1, s).second && (gr.get_vertex_info(s-1).rpos==gr.get_vertex_info(s).lpos))
        {
            k = gx.locate_vertex(vi.lpos-1);
        }
        if(k < 0 || k == gx.num_vertices()-1 || gx.edge(k, gx.num_vertices()-1).second) continue;
        double new_loss = 0;
        if(gx.edge(k, k+1).second && (gx.get_vertex_info(k).rpos == gx.get_vertex_info(k+1).lpos))
        {
            new_loss = gx.get_out_weights(k)-gx.get_edge_weight(gx.edge(k, k+1).first);
        }
        else
        {
            new_loss = gx.get_out_weights(k);
        }
        if(cfg.verbose >= 2) printf("End vertex %d(gx=%d, vertex=%d) boundary_loss = %.2lf\n", s, sample_id, k, new_loss);
        if(sample_id == -1 && pos_type == 1)
            vi.boundary_merged_loss += new_loss; 
        else
        {
            if(pos_type == 1)
                vi.boundary_loss1 += new_loss; 
            else if(pos_type == 2)
                vi.boundary_loss2 += new_loss; 
            else if(pos_type == 3)
                vi.boundary_loss3 += new_loss; 
        }
        //gr.set_vertex_info(s, vi);

    }

    //bridge bonus
    /*pei = gr.out_edges(0);
    for(it = pei.first; it != pei.second; it++)
    {
        edge_descriptor e = *it;
        int s = e->source();
        int t = e->target();
        assert(s == 0);

        if(t-1 <= 0) continue;
        if(!gr.edge(t-1, gr.num_vertices()-1).second) continue;
        vertex_info vi1 = gr.get_vertex_info(t-1);
        int k1 = gx.locate_vertex(vi1.rpos -1);
        vertex_info vi2 = gr.get_vertex_info(t);
        int k2 = gx.locate_vertex(vi2.lpos);
        if(k1 > 0 && k2 > 0 && k1 == k2)//possible bridge
        {
            if(gx.edge(0, k1).second || gx.edge(k1, gx.num_vertices()-1).second) continue;
            vi1.bridge_bonus += gx.get_vertex_weight(k1);
            vi2.bridge_bonus += gx.get_vertex_weight(k2);
            gr.set_vertex_info(t-1, vi1);
            gr.set_vertex_info(t, vi2);
            if(cfg.verbose >= 2) printf("Bridge %d-%d(gx:%d) bonus = %.2lf, %.2lf\n", t-1, t, k1, vi1.bridge_bonus, vi2.bridge_bonus);

        }
    }*/

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
	/*
	for(int k = 0; k < cfg.assembly_repeats; k++)
	{
		boost::asio::post(pool, [gx, hx, k, sid, pa, &mt, &tsp, &tm] {
		*/

        //splice_graph gr(gx);
        hyper_set hs(hx);

        transcript_set ts(gx.chrm, tm.rid, pa.min_single_exon_clustering_overlap);

        //printf("A: tm.rid = %d, ts.rid = %d, this->rid = %d\n", tm.rid, ts.rid, this->rid);

		int k = 0;
        gx.gid = gx.gid + "." + tostring(k);
        scallop sx(gx, hs, pa, k == 0 ? false : true);
        sx.assemble();

        int z = 0;
        for(int i = 0; i < sx.trsts.size(); i++)
        {
            transcript &t = sx.trsts[i];
            z++;
            t.RPKM = 0;
            ts.add(t, 1, sid, TRANSCRIPT_COUNT_ADD_COVERAGE_ADD);
        }

        if(pa.verbose >= 2) printf("assemble %s: %d transcripts, graph with %lu vertices and %lu edges\n", gx.gid.c_str(), z, gx.num_vertices(), gx.num_edges());
        if(gx.num_vertices() >= 1000) printf("assemble %s: %d transcripts, large graph with %lu vertices and %lu edges\n", gx.gid.c_str(), z, gx.num_vertices(), gx.num_edges());

        mt.lock();
        //tsp.tsets.push_back(ts);
        //printf("B: tm.rid = %d, ts.rid = %d, this->rid = %d\n", tm.rid, ts.rid, this->rid);
        tm.add(ts, TRANSCRIPT_COUNT_ADD_COVERAGE_ADD);
        mt.unlock();

		/*
		});
	}
	*/

	return 0;
}
