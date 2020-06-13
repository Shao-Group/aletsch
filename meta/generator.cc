/*
Part of meta-scallop
(c) 2020 by Mingfu Shao, The Pennsylvania State University
See LICENSE for licensing.
*/

#include <cstdio>
#include <cassert>
#include <sstream>
#include "boost/pending/disjoint_sets.hpp"
#include <boost/asio/post.hpp>
#include <boost/asio/thread_pool.hpp>
#include <boost/pending/disjoint_sets.hpp>

#include "constants.h"
#include "parameters.h"
#include "bundle.h"
#include "generator.h"
#include "previewer.h"
#include "bridge_solver.h"
#include "graph_builder.h"
#include "graph_cluster.h"
#include "graph_reviser.h"
#include "essential.h"
#include "hyper_set.h"
#include "assembler.h"

generator::generator(sample_profile &s, vector<combined_graph> &v, vector<transcript> &t, const parameters &c, int tid)
	: vcb(v), trsts(t), cfg(c), sp(s), target_id(tid)
{
	index = 0;
	qlen = 0;
	qcnt = 0;
	sp.open_align_file();
}

generator::~generator()
{
	sp.close_align_file();
}

int generator::resolve()
{
	if(target_id < 0) return 0;

	int index = 0;
	bundle bb1;
	bundle bb2;

	int hid = 0;
    bam1_t *b1t = bam_init1();

	hts_itr_t *iter = sp.iters[target_id];
	if(iter == NULL) return 0;

	while(sam_itr_next(sp.sfn, iter, b1t) >= 0)
	{
		bam1_core_t &p = b1t->core;

		if((p.flag & 0x4) >= 1) continue;											// read is not mapped
		if((p.flag & 0x100) >= 1 && cfg.use_second_alignment == false) continue;	// secondary alignment
		if(p.n_cigar > cfg.max_num_cigar) continue;									// ignore hits with more than max-num-cigar types
		if(p.qual < cfg.min_mapping_quality) continue;								// ignore hits with small quality
		if(p.n_cigar < 1) continue;													// should never happen

		hit ht(b1t, hid++);
		ht.set_splices(b1t);
		ht.set_tags(b1t);
		ht.set_strand(sp.library_type);

		//ht.print();

		qlen += ht.qlen;
		qcnt += 1;

		// truncate
		if(bb1.hits.size() >= 1 && (ht.tid != bb1.tid || ht.pos > bb1.rpos + cfg.min_bundle_gap))
		{
			generate(bb1, index);
			bb1.clear();
			index++;
		}

		if(bb2.hits.size() >= 1 && (ht.tid != bb2.tid || ht.pos > bb2.rpos + cfg.min_bundle_gap))
		{
			generate(bb2, index);
			bb2.clear();
			index++;
		}

		// add hit
		if(cfg.uniquely_mapped_only == true && ht.nh != 1) continue;
		if(sp.library_type != UNSTRANDED && ht.strand == '+' && ht.xs == '-') continue;
		if(sp.library_type != UNSTRANDED && ht.strand == '-' && ht.xs == '+') continue;
		if(sp.library_type != UNSTRANDED && ht.strand == '.' && ht.xs != '.') ht.strand = ht.xs;
		if(sp.library_type != UNSTRANDED && ht.strand == '+') bb1.add_hit_intervals(ht, b1t);
		if(sp.library_type != UNSTRANDED && ht.strand == '-') bb2.add_hit_intervals(ht, b1t);
		if(sp.library_type == UNSTRANDED) bb1.add_hit_intervals(ht, b1t);
		/*
		if(sp.library_type == UNSTRANDED && ht.xs == '+') bb1.add_hit_intervals(ht, b1t);
		if(sp.library_type == UNSTRANDED && ht.xs == '-') bb2.add_hit_intervals(ht, b1t);
		if(sp.library_type == UNSTRANDED && ht.xs == '.') bb1.add_hit_intervals(ht, b1t);
		if(sp.library_type == UNSTRANDED && ht.xs == '.') bb2.add_hit_intervals(ht, b1t);
		*/
	}

    bam_destroy1(b1t);

	generate(bb1, index++);
	generate(bb2, index++);
	bb1.clear();
	bb2.clear();

	return 0;
}

int generator::generate(bundle &bb, int index)
{
	if(bb.tid < 0) return 0;

	bool store_hits = true;
	if(cfg.output_bridged_bam_dir == "") store_hits = false;
	if(sp.data_type != PAIRED_END) store_hits = false;

	if(bb.hits.size() < cfg.min_num_hits_in_bundle)
	{
		if(store_hits == false || bb.hits.size() <= 0) return 0;

		sp.open_bridged_bam(cfg.output_bridged_bam_dir);
		for(int i = 0; i < bb.hits.size(); i++)
		{
			hit &h = bb.hits[i];
			bam1_t b1t;
			bool b = build_bam1_t(b1t, h, bb.hits[i].spos);
			if(b == true) bam_write1(sp.bridged_bam, &(b1t));
			assert(b1t.data != NULL);
			delete b1t.data;
		}
		sp.close_bridged_bam();
		return 0;
	}

	char buf[1024];
	strcpy(buf, sp.hdr->target_name[bb.tid]);
	bb.chrm = string(buf);
	//bb.compute_strand(sp.library_type);

	splice_graph gr;
	graph_builder gb(bb, cfg);
	gb.build(gr);
	gr.extend_strands();
	gr.build_vertex_index();

	if(sp.data_type == PAIRED_END || sp.data_type == SINGLE_END)
	{
		revise_splice_graph_full(gr, cfg);
	}

	/*
	printf("-----------------------------\n");
	gr.print();
	printf("\n");
	*/

	phase_set ps;
	vector<pereads_cluster> ub;

	if(sp.data_type != PAIRED_END)
	{
		vector<bool> paired(bb.hits.size(), false);
		build_phase_set_from_unpaired_reads(ps, gr, bb.hits, paired);
	}
	else
	{
		vector<pereads_cluster> vc;
		graph_cluster gc(gr, bb.hits, cfg.max_reads_partition_gap, store_hits);
		gc.build_pereads_clusters(vc);

		vector<bool> paired = gc.get_paired();
		build_phase_set_from_unpaired_reads(ps, gr, bb.hits, paired);

		bridge_solver bs(gr, vc, cfg, sp.insertsize_low, sp.insertsize_high);
		bs.build_phase_set(ps);

		bs.collect_unbridged_clusters(ub); 
		if(store_hits == true)
		{
			sp.open_bridged_bam(cfg.output_bridged_bam_dir);
			write_unpaired_reads(sp.bridged_bam, bb.hits, paired);
			assert(vc.size() == bs.opt.size());
			for(int k = 0; k < vc.size(); k++)
			{
				if(bs.opt[k].type >= 0) write_bridged_pereads_cluster(sp.bridged_bam, vc[k], bs.opt[k].whole);
				else write_unbridged_pereads_cluster(sp.bridged_bam, vc[k]);
			}
			sp.close_bridged_bam();
		}
		for(int i = 0; i < vc.size(); i++) vc[i].clear();
	}

	//printf("single-bridge, combined = %d, ", 1); bs.print();

	vector<splice_graph> grv;
	vector<phase_set> hsv;
	vector< vector<pereads_cluster> > ubv;
	partition(gr, ps, ub, grv, hsv, ubv);

	assert(grv.size() == hsv.size());
	assert(grv.size() == ubv.size());

	for(int k = 0; k < grv.size(); k++)
	{
		string gid = "gene." + tostring(sp.sample_id) + "." + tostring(index) + "." + tostring(k);
		grv[k].gid = gid;

		bool b = regional(grv[k], hsv[k], ubv[k]);
		if(b == true) continue;

		/*
		b = assemble(grv[k], hsv[k], ubv[k], tlock);
		if(b == true) continue;
		*/

		combined_graph cb(cfg);
		cb.sid = sp.sample_id;
		cb.gid = gid;
		cb.build(grv[k], hsv[k], ubv[k]);
		//cb.refine_junctions();

		// print
		/*
		grv[k].print();
		cb.print(k);
		printf("\n");
		*/

		vcb.push_back(cb);
	}

	return 0;
}

bool generator::regional(splice_graph &gr, phase_set &ps, vector<pereads_cluster> &vc)
{
	bool all_regional = true;
	for(int i = 1; i < gr.num_vertices() - 1; i++)
	{
		if(gr.get_vertex_info(i).regional == false) all_regional = false;
		if(all_regional == false) break;
	}

	if(all_regional == false && gr.num_edges() >= 1) return false;

	if(cfg.output_bridged_bam_dir != "" && vc.size() >= 1)
	{
		sp.open_bridged_bam(cfg.output_bridged_bam_dir);
		for(int k = 0; k < vc.size(); k++)
		{
			write_unbridged_pereads_cluster(sp.bridged_bam, vc[k]);
			vc[k].clear();
		}
		sp.close_bridged_bam();
	}

	//gr.print(); printf("above graph is a regional graph\n\n");

	return true;
}

bool generator::assemble(splice_graph &gr, phase_set &ps, vector<pereads_cluster> &vc, mutex &tlock)
{
	//if(gr.num_vertices() >= 3 && gr.num_vertices() <= 1000) return false;
	if(gr.num_vertices() <= 1000) return false;

	vector<transcript> vt;
	assembler asmb(cfg);
	asmb.assemble(gr, ps, vt);

	tlock.lock();
	trsts.insert(trsts.end(), vt.begin(), vt.end());
	tlock.unlock();

	if(cfg.output_bridged_bam_dir != "" && vc.size() >= 1)
	{
		sp.open_bridged_bam(cfg.output_bridged_bam_dir);
		for(int k = 0; k < vc.size(); k++)
		{
			write_unbridged_pereads_cluster(sp.bridged_bam, vc[k]);
			vc[k].clear();
		}
		sp.close_bridged_bam();
	}

	return true;
}

int generator::partition(splice_graph &gr, phase_set &ps, vector<pereads_cluster> &ub, vector<splice_graph> &grv, vector<phase_set> &psv, vector< vector<pereads_cluster> > &ubv)
{
	int n = gr.num_vertices();

	vector<int> rank(n, -1);
	vector<int> parent(n, -1);

	disjoint_sets<int*, int*> ds(&rank[0], &parent[0]);
	for(int k = 0; k < n; k++) ds.make_set(k);

	// group with edges in gr
	PEEI pei = gr.edges();
	for(edge_iterator it = pei.first; it != pei.second; it++)
	{
		edge_descriptor e = (*it);
		int s = e->source();
		int t = e->target();
		if(s == 0) continue;
		if(t == n - 1) continue;
		ds.union_set(s, t);
	}

	// group with phase_set, do not use
	/*
	for(MVII::const_iterator it = ps.pmap.begin(); it != ps.pmap.end(); it++)
	{
		const vector<int> &v = it->first;
		assert(v.size() % == 0);
		assert(v.size() >= 2);
		
		assert(gr.lindex.find(v.front()) != gr.lindex.end());
		assert(gr.rindex.find(v.back()) != gr.rindex.end());
		int kl = gr.lindex[v.front()];
		int kr = gr.rindex[v.back()];
		ds.union_set(kl, kr);
	}
	*/

	// group with unbridged pairs
	for(int i = 0; i < ub.size(); i++)
	{
		/*
		int32_t p1 = ub[i].bounds[1] - 1;
		int32_t p2 = ub[i].bounds[2] - 0;
		int x = gr.locate_vertex(p1);
		int y = gr.locate_vertex(p2);
		*/
		assert(gr.rindex.find(ub[i].extend[1]) != gr.rindex.end());
		assert(gr.lindex.find(ub[i].extend[2]) != gr.lindex.end());
		int x = gr.rindex[ub[i].extend[1]];
		int y = gr.lindex[ub[i].extend[2]];
		ds.union_set(x, y);
	}

	// create connected components
	vector< set<int> > vv;
	map<int, int> m;
	for(int i = 1; i < n - 1; i++)
	{
		int p = ds.find_set(i);
		if(m.find(p) == m.end())
		{
			m.insert(pair<int, int>(p, vv.size()));
			set<int> v;
			v.insert(i);
			vv.push_back(v);
		}
		else
		{
			int k = m[p];
			vv[k].insert(i);
		}
	}

	// print
	/*
	for(int i = 0; i < vv.size(); i++)
	{
		printf("component %d: ", i);
		vector<int> z(vv[i].begin(), vv[i].end());
		printv(z);
		printf("\n");
	}
	*/

	vector< map<int, int> > vm;
	for(int k = 0; k < vv.size(); k++)
	{
		map<int, int> a2b;
		transform_vertex_set_map(vv[k], a2b);
		vm.push_back(a2b);
	}

	grv.resize(vv.size());
	for(int k = 0; k < vv.size(); k++)
	{
		build_child_splice_graph(gr, grv[k], vm[k]);
	}

	psv.resize(vv.size());
	for(MVII::const_iterator it = ps.pmap.begin(); it != ps.pmap.end(); it++)
	{
		const vector<int> &v = it->first;
		assert(v.size() % 2 == 0);
		if(v.size() <= 1) continue;
		
		assert(gr.lindex.find(v.front()) != gr.lindex.end());
		assert(gr.rindex.find(v.back()) != gr.rindex.end());
		int j = gr.lindex[v.front()];
		int p = ds.find_set(j);
		assert(m.find(p) != m.end());
		int k = m[p];
		assert(k >= 0 && k < vv.size());
		psv[k].add(v, it->second);
	}

	ubv.resize(vv.size());
	for(int i = 0; i < ub.size(); i++)
	{
		int x = gr.rindex[ub[i].extend[1]];
		int p = ds.find_set(x);
		int k = m[p];
		assert(k >= 0 && k < vv.size());
		ubv[k].push_back(std::move(ub[i]));
	}

	return 0;
}
