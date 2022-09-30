/*
Part of aletsch
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
#include "essential.h"
#include "hyper_set.h"
#include "assembler.h"

generator::generator(sample_profile &s, vector<bundle> &v, transcript_set &t, const parameters &c, int tid, int rid)
	: vcb(v), ts(t), cfg(c), sp(s), target_id(tid), region_id(rid)
{
	index = 0;
	sp.open_align_file();
}

generator::~generator()
{
	sp.close_align_file();
}

int generator::resolve()
{
	if(target_id < 0 || region_id < 0) return 0;

	int index = 0;
	bundle_base bb1;
	bundle_base bb2;

	int hid = 0;
    bam1_t *b1t = bam_init1();

	hts_itr_t *iter = sp.iters[target_id][region_id];
	if(iter == NULL) return 0;

	int32_t start1 = sp.start1[target_id][region_id];
	int32_t start2 = sp.start2[target_id][region_id];
	int32_t new_start1 = start1;
	int32_t new_start2 = start2;

	bool term1 = false, term2 = false;
	while(sam_itr_next(sp.sfn, iter, b1t) >= 0)
	{
		bam1_core_t &p = b1t->core;

		if(p.pos < sp.region_partition_length * region_id) continue;
		if(p.pos < start1 && p.pos < start2) continue;

		if((p.flag & 0x4) >= 1) continue;											// read is not mapped
		if((p.flag & 0x100) >= 1 && cfg.use_second_alignment == false) continue;	// secondary alignment
		if(p.n_cigar > cfg.max_num_cigar) continue;									// ignore hits with more than max-num-cigar types
		if(p.qual < cfg.min_mapping_quality) continue;								// ignore hits with small quality
		if(p.n_cigar < 1) continue;													// should never happen

		hit ht(b1t, hid++);
		ht.set_tags(b1t);
		ht.set_strand(sp.library_type);
		//ht.print();

		// truncate
		if(bb1.hits.size() >= 1 && (ht.tid != bb1.tid || ht.pos > bb1.rpos + cfg.min_bundle_gap))
		{
			generate(bb1, index);
			bb1.clear();
			index++;
			if(ht.pos >= sp.region_partition_length * (1 + region_id)) term1 = true;
		}
		if(bb1.hits.size() <= 0 && ht.pos >= sp.region_partition_length * (1 + region_id)) term1 = true;
		if(term1 == true) new_start1 = ht.pos;

		if(bb2.hits.size() >= 1 && (ht.tid != bb2.tid || ht.pos > bb2.rpos + cfg.min_bundle_gap))
		{
			generate(bb2, index);
			bb2.clear();
			index++;
			if(ht.pos >= sp.region_partition_length * (1 + region_id)) term2 = true;
		}
		if(bb2.hits.size() <= 0 && ht.pos >= sp.region_partition_length * (1 + region_id)) term2 = true;
		if(term2 == true) new_start2 = ht.pos;

		if(term1 == true && term2 == true) break;

		// add hit
		if(cfg.uniquely_mapped_only == true && ht.nh != 1) continue;
		if(sp.library_type != UNSTRANDED && ht.strand == '+' && ht.xs == '-') continue;
		if(sp.library_type != UNSTRANDED && ht.strand == '-' && ht.xs == '+') continue;
		//if(sp.library_type == UNSTRANDED && sp.bam_with_xs == 1 && ht.xs == '.') continue;
		if(sp.library_type != UNSTRANDED && ht.strand == '.' && ht.xs != '.') ht.strand = ht.xs;

		if(sp.library_type != UNSTRANDED && ht.strand == '+' && ht.pos >= start1 && term1 == false) bb1.add_hit_intervals(ht, b1t);
		if(sp.library_type != UNSTRANDED && ht.strand == '-' && ht.pos >= start2 && term2 == false) bb2.add_hit_intervals(ht, b1t);

		// TODO, handle unstranded case
		assert(sp.library_type != UNSTRANDED);
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

	if(region_id < sp.start1[target_id].size() - 1) sp.start1[target_id][region_id + 1] = new_start1;
	if(region_id < sp.start2[target_id].size() - 1) sp.start2[target_id][region_id + 1] = new_start2;

	return 0;
}

int generator::generate(bundle_base &bb, int index)
{
	if(bb.tid < 0) return 0;
	char buf[1024];
	strcpy(buf, sp.hdr->target_name[bb.tid]);

	bundle bd(cfg, sp, std::move(bb));
	bd.chrm = string(buf);
	bd.gid = "gene." + tostring(sp.sample_id) + "." + tostring(index);
	bd.build_fragments();
	bd.bridge();
	//bd.filter_multialigned_hits();

	// TODO, storing reads
	// TODO, don't keep bridged reads

	vcb.push_back(std::move(bd));
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

bool generator::assemble_single(splice_graph &gr, phase_set &ps, vector<pereads_cluster> &vc)
{
	transcript t;
	bool b = build_single_exon_transcript(gr, t);
	if(b == false) return false;

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

	//if(t.coverage < cfg.min_single_exon_transcript_coverage) return true;
	if(t.length() < cfg.min_single_exon_transcript_length) return true;

	ts.add(t, 2, sp.sample_id, TRANSCRIPT_COUNT_ADD_COVERAGE_ADD);
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
	for(MVII::iterator it = ps.pmap.begin(); it != ps.pmap.end(); it++)
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
