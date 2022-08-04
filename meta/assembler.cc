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
		// print
		printf("\n");
		for(int k = 0; k < gv.size(); k++) gv[k]->print(k);

		vector<vector<PID>> sim;
		build_similarity(gv, sim);
		bridge_pairwise(gv, sim);
		//refine_pairwise(gv, sim);
		assemble(gv, ts, instance);
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

		fix_missing_edges(gr, gx);

		phase_set ps;
		bd.build_phase_set(ps, gr);
		px.combine(ps);

		assemble(gr, ps, ts, bd.sp.sample_id);
	}

	// assemble combined instance
	assemble(gx, px, ts, -1);
	return 0;
}

int assembler::refine_pairwise(vector<bundle*> gv, vector<vector<PID>> &sim)
{
	for(int i = 0; i < gv.size(); i++)
	{
		for(int k = 0; k < sim[i].size(); k++)
		{
			int j = sim[i][k].first;
			vector<bundle*> v;
			v.push_back(gv[i]);
			v.push_back(gv[j]);
			printf("trying to refine pair %d and %d, sim = %.2lf, samples %d and %d, files %s and %s\n", 
					i, j, sim[i][k].second, gv[i]->sp.sample_id, gv[j]->sp.sample_id, gv[i]->sp.align_file.c_str(), gv[j]->sp.align_file.c_str());
			refine(v);
		}
	}
	for(int i = 0; i < gv.size(); i++)
	{
		printf("print borrowed paths for bundle %d\n", i);
		gv[i]->print_borrowed_paths();
		gv[i]->digest_borrowed_paths();
	}
	return 0;
}

int assembler::refine(vector<bundle*> gv)
{
	assert(gv.size() >= 2);

	// construct combined bundle
	bundle cb(cfg, gv[0]->sp);
	cb.copy_meta_information(*(gv[0]));
	for(int k = 0; k < gv.size(); k++) cb.combine(*(gv[k]));

	// construct combined graph
	splice_graph gr;
	transform(cb, gr, false);

	// refine each individual bundle
	for(int k = 0; k < gv.size(); k++)
	{
		refine(gv[k], gr);
	}
	return 0;
}

/*
int assembler::refine_pairwise(bundle &cx, bundle &cy)
{
	// combined bundle
	bundle cb(cfg, cx.sp);
	cb.copy_meta_information(cx);
	cb.combine(cx);
	cb.combine(cy);
	cb.set_gid(0, 0);		// TODO, proper gid

	// combined splice graph
	splice_graph gr;
	transform(cb, gr, false);

	// individual graphs
	splice_graph gx;
	splice_graph gy;
	transform(cx, gx, false);
	transform(cy, gy, false);

	refine_pairwise(gx, gr);
	refine_pairwise(gy, gr);
	return 0;
}
*/

int assembler::refine(bundle *bd, splice_graph &gr)
{
	splice_graph gx;
	transform(*bd, gx, false);

	printf("refine gr and gx: strand = %c and %c\n", gx.strand, gr.strand);
	if(gx.strand != gr.strand) return 0;
	int strand = 0;
	if(gr.strand == '+') strand = 1;
	else if(gr.strand == '-') strand = 2;
	if(strand == 0) return 0; // TODO

	// build index for all starting and ending positions
	int m = gx.num_vertices() - 1;
	set<int32_t> sset;
	set<int32_t> tset;
	PEEI pei = gx.out_edges(0);
	for(edge_iterator it = pei.first; it != pei.second; it++)
	{
		int s = (*it)->source();
		int t = (*it)->target();
		assert(s == 0);
		int32_t z = gx.get_vertex_info(t).lpos;
		sset.insert(z);
	}
	pei = gx.in_edges(m);
	for(edge_iterator it = pei.first; it != pei.second; it++)
	{
		int s = (*it)->source();
		int t = (*it)->target();
		assert(t == m);
		int32_t z = gx.get_vertex_info(s).rpos;
		tset.insert(z);
	}

	// DP starting from the source 0
	int n = gr.num_vertices() - 1;
	vector<pereads_cluster> vc;
	bridge_solver bs0(gr, vc, cfg);
	vector<vector<entry>> table0; 
	bs0.dynamic_programming(0, n, table0, strand);

	// check all starting vertices of gx
	pei = gx.out_edges(0);
	for(edge_iterator it = pei.first; it != pei.second; it++)
	{
		int s = (*it)->source();
		int t = (*it)->target();
		assert(s == 0);

		// TODO, do not refine isolated vertices
		if(gx.edge(t, gx.num_vertices() - 1).second == true) continue;

		int32_t z = gx.get_vertex_info(t).lpos;
		int v = gr.locate_vertex(z);

		if(v < 0 || v >= gr.num_vertices()) continue;

		// if v is also a starting vertex continue;
		//printf("start vertex of gx = %d, locating to gr.v = %d, start = %c, end = %c\n", t, v, gr.edge(0, v).second ? 'T' : 'F', gr.edge(v, gr.num_vertices() - 1).second ? 'T' : 'F');

		if(gr.edge(0, v).second == true) continue;
		if(gr.edge(v, gr.num_vertices() - 1).second == true) continue;

		vector<vector<int>> pb = bs0.trace_back(v, table0);

		for(int j = 0; j < pb.size(); j++)
		{
			bridge_path p;
			p.score = table0[v][j].stack.front();
			p.stack = table0[v][j].stack;
			p.v = pb[j];
			build_intron_coordinates_from_path(gr, p.v, p.chain);

			printf("bridging chain %d -> %d: ", gr.get_vertex_info(0).lpos, z);
			printv(p.chain);
			printf("\n");
			printf("chain.score = %.2lf, chain.stack = ", p.score);
			printv(p.stack);
			printf("\n");

			// annotate path
			vector<int32_t> pt;
			pt.push_back(gr.get_vertex_info(0).lpos);
			pt.insert(pt.end(), p.chain.begin(), p.chain.end());
			pt.push_back(z);

			if(pt[0] == pt[1])
			{
				pt.erase(pt.begin());
				pt.erase(pt.begin());
			}
			vector<int32_t> vv, nn, pp;
			annotate_path(gx, pt, vv, nn, pp);
			assert(vv.size() == nn.size());
			assert(vv.size() % 2 == 0);
			for(int k = 0; k < vv.size() / 2; k++)
			{
				printf(" region: %9d - %9d, length = %5d, type = %1d, annotate = %2d, ltype = %d, rtype = %d, lbound = %c/%c, rbound = %c/%c\n", 
						vv[k * 2 + 0], vv[k * 2 + 1], vv[k * 2 + 1] - vv[k * 2 + 0], nn[k * 2 + 0], nn[k * 2 + 1], pp[k * 2 + 0], pp[k * 2 + 1],
						sset.find(vv[k*2+0]) == sset.end() ? 'F' : 'T', tset.find(vv[k*2+0]) == tset.end() ? 'F' : 'T', sset.find(vv[k*2+1]) == sset.end() ? 'F' : 'T', tset.find(vv[k*2+1]) == tset.end() ? 'F' : 'T');
			}
			printf("\n");

			//p.chain = filter_pseudo_introns(p.chain);
			//piers[b].bridges.push_back(p);

			//determine if it is correct

			int m = vv.size() / 2;
			bool accept = true;
			vector<int32_t> bpath;
			double bweight = -1;
			
			if(m >= 3 && nn[m*2-1] == -1 && nn[m*2-2] == 1 && nn[m*2-3] == -1 && nn[m*2-4] == 2 && nn[m*2-5] == -1 && nn[m*2-6] == 1)
			{
				if(m == 3) accept = false;
				for(int k = m - 3; k >= 2; k--)
				{
					if(nn[k*2-1] == -1) accept = false;
					if(accept == false) break;
				}
				int32_t len1 = vv[m*2-1] - vv[m*2-2];
				int32_t len2 = vv[m*2-3] - vv[m*2-4];
				int32_t len3 = vv[m*2-5] - vv[m*2-6];

				if(len2 < len1 || len2 < len3) accept = false;
				if(len1 > 100 || len3 > 100) accept = false;

				if(accept == true)
				{
					printf("ACCEPT PATH, type = (-1,-1,-1), (1, 2, 1)\n");
					bpath.push_back(0 + vv[m*2-6]);
					bpath.push_back(0 + vv[m*2-5]);
					bpath.push_back(0 - vv[m*2-4]);
					bpath.push_back(0 - vv[m*2-3]);
					bpath.push_back(0 + vv[m*2-2]);
					bpath.push_back(0 + vv[m*2-1]);
					bweight = p.score;
				}
				else printf("REJECT PATH, type = (-1,-1,-1), (1, 2, 1)\n");
			}
			else if(m >= 3 && nn[m*2-1] == 1 && nn[m*2-2] == 1 && nn[m*2-3] == -1 && nn[m*2-4] == 2 && nn[m*2-5] == -1 && nn[m*2-6] == 1)
			{
				for(int k = m - 3; k >= 2; k--)
				{
					if(nn[k*2-1] == -1) accept = false;
					if(accept == false) break;
				}
				int32_t len2 = vv[m*2-3] - vv[m*2-4];
				int32_t len3 = vv[m*2-5] - vv[m*2-6];

				if(len2 < len3) accept = false;
				if(len3 > 100) accept = false;

				if(accept == true)
				{
					printf("ACCEPT PATH, type = (1, -1, -1), (1, 2, 1)\n");
					bpath.push_back(0 - vv[m*2-4]);
					bpath.push_back(0 - vv[m*2-3]);
					bpath.push_back(0 + vv[m*2-2]);
					bpath.push_back(0 + vv[m*2-1]);
					bweight = p.score;
				}
				else printf("REJECT PATH, type = (1, -1, -1), (1, 2, 1)\n");
			}
			else if(m >= 2 && nn[m*2-1] == -1 && nn[m*2-2] == 2 && nn[m*2-3] == -1 && nn[m*2-4] == 1)
			{
				for(int k = m - 2; k >= 2; k--)
				{
					if(nn[k*2-1] == -1) accept = false;
					if(accept == false) break;
				}
				int32_t len1 = vv[m*2-1] - vv[m*2-2];
				int32_t len2 = vv[m*2-3] - vv[m*2-4];

				if(len1 < len2) accept = false;
				if(len2 > 100) accept = false;

				if(accept == true)
				{
					bpath.push_back(0 + vv[m*2-4]);
					bpath.push_back(0 + vv[m*2-3]);
					bpath.push_back(0 - vv[m*2-2]);
					bpath.push_back(0 - vv[m*2-1]);
					bweight = p.score;
					printf("ACCEPT PATH, type = (-1, -1), (2, 1)\n");
				}
				else printf("REJECT PATH, type = (-1, -1), (2, 1)\n");
			}
			else if(m >= 2 && nn[m*2-1] == -1 && nn[m*2-2] == 1 && nn[m*2-3] == -1 && nn[m*2-4] == 2)
			{
				for(int k = m - 2; k >= 2; k--)
				{
					if(nn[k*2-1] == -1) accept = false;
					if(accept == false) break;
				}
				int32_t len1 = vv[m*2-1] - vv[m*2-2];
				int32_t len2 = vv[m*2-3] - vv[m*2-4];

				if(len2 < len1) accept = false;
				if(len1 > 100) accept = false;

				if(accept == true)
				{
					bpath.push_back(0 - vv[m*2-4]);
					bpath.push_back(0 - vv[m*2-3]);
					bpath.push_back(0 + vv[m*2-2]);
					bpath.push_back(0 + vv[m*2-1]);
					bweight = p.score;
					printf("ACCEPT PATH, type = (-1, -1), (1, 2)\n");
				}
				else printf("REJECT PATH, type = (-1, -1), (1, 2)\n");
			}
			else if(m >= 1 && nn[m*2-1] == -1 && nn[m*2-2] == 2)
			{
				for(int k = m - 1; k >= 2; k--)
				{
					if(nn[k*2-1] == -1) accept = false;
					if(accept == false) break;
				}
				if(accept == true) 
				{
					printf("ACCEPT PATH, type = (-1), (2)\n");
					bpath.push_back(0 - vv[m*2-2]);
					bpath.push_back(0 - vv[m*2-1]);
					bweight = p.score;
				}
				else printf("REJECT PATH, type = (-1), (2)\n");
			}
			else
			{
				accept = false;
				printf("REJECT PATH, other type\n");
			}

			if(accept == true)
			{
				bd->save_borrowed_path(bpath, bweight);
			}

			break;
		}
	}
	
	return 0;
}

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
		printf("fixing starting boundary t = %d-%d using u = %d-%d, v = %d-%d, gap = %d, wt = %.1lf, wuv = %.1lf\n", 
				vt.lpos, vt.rpos, vu.lpos, vu.rpos, vv.lpos, vv.rpos, gap, wt, wuv);
	}

	return 0;
}

int assembler::bridge_pairwise(vector<bundle*> gv, vector<vector<PID>> &sim)
{
	for(int i = 0; i < gv.size(); i++)
	{
		int c = gv[i]->count_unbridged();
		if(c <= 0) continue;

		for(int k = 0; k < sim[i].size(); k++)
		{
			int j = sim[i][k].first;
			vector<bundle*> v;
			v.push_back(gv[i]);
			v.push_back(gv[j]);
			bridge(v);

			//printf("trying to bridge pair %d and %d, sim = %.2lf, samples %d and %d, files %s and %s\n", i, j, sim[i][k].second, gv[i]->sp.sample_id, gv[j]->sp.sample_id, gv[i]->sp.align_file.c_str(), gv[j]->sp.align_file.c_str());
		}
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
			cnt2 += bd.update_bridges(vc[j].frlist, bs.opt[j].chain);
			//vc[j].print(j);
		}

		//if(cfg.verbose >= 2) 
		//printf("gid %s: further bridge %d / %lu clusters, %d / %d fragments\n", bd.gid.c_str(), cnt1, vc.size(), cnt2, unbridged);
	}
	return 0;
}

int assembler::assemble(splice_graph &gx, phase_set &px, transcript_set &ts, int sid)
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
