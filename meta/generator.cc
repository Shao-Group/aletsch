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

generator::generator(sample_profile &s, vector<combined_graph> &v, const parameters &c)
	: vcb(v), cfg(c), sp(s)
{
	previewer pre(cfg, sp);
	pre.infer_library_type();
	pre.infer_insertsize();

    sfn = sam_open(sp.file_name.c_str(), "r");
    hdr = sam_hdr_read(sfn);
    b1t = bam_init1();
	index = 0;
	qlen = 0;
	qcnt = 0;
}

generator::~generator()
{
    bam_destroy1(b1t);
    bam_hdr_destroy(hdr);
    sam_close(sfn);
}

int generator::resolve()
{
	int num_threads = 0;
	if(cfg.single_sample_multiple_threading) num_threads = 1;

	boost::asio::thread_pool pool(num_threads);				// thread pool
	mutex mylock;											// lock for trsts

	int index = 0;
	bundle *bb1 = new bundle();
	bundle *bb2 = new bundle();

	int hid = 0;

    while(sam_read1(sfn, hdr, b1t) >= 0)
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

		// TODO for test
		//if(ht.tid >= 1) break;
		//if(index >= 100) break;

		qlen += ht.qlen;
		qcnt += 1;

		// truncate
		if(bb1->hits.size() >= 1 && (ht.tid != bb1->tid || ht.pos > bb1->rpos + cfg.min_bundle_gap))
		{
			if(cfg.single_sample_multiple_threading == false) this->generate(bb1, mylock, index);
			else boost::asio::post(pool, [this, &mylock, bb1, index]{ this->generate(bb1, mylock, index); });
			bb1 = new bundle();
			index++;
		}

		if(bb2->hits.size() >= 1 && (ht.tid != bb2->tid || ht.pos > bb2->rpos + cfg.min_bundle_gap))
		{
			if(cfg.single_sample_multiple_threading == false) this->generate(bb2, mylock, index);
			else boost::asio::post(pool, [this, &mylock, bb2, index]{ this->generate(bb2, mylock, index); });
			bb2 = new bundle();
			index++;
		}

		// add hit
		if(cfg.uniquely_mapped_only == true && ht.nh != 1) continue;
		if(sp.library_type != UNSTRANDED && ht.strand == '+' && ht.xs == '-') continue;
		if(sp.library_type != UNSTRANDED && ht.strand == '-' && ht.xs == '+') continue;
		if(sp.library_type != UNSTRANDED && ht.strand == '.' && ht.xs != '.') ht.strand = ht.xs;
		if(sp.library_type != UNSTRANDED && ht.strand == '+') bb1->add_hit_intervals(ht, b1t);
		if(sp.library_type != UNSTRANDED && ht.strand == '-') bb2->add_hit_intervals(ht, b1t);
		if(sp.library_type == UNSTRANDED) bb1->add_hit_intervals(ht, b1t);
		/*
		if(sp.library_type == UNSTRANDED && ht.xs == '+') bb1->add_hit_intervals(ht, b1t);
		if(sp.library_type == UNSTRANDED && ht.xs == '-') bb2->add_hit_intervals(ht, b1t);
		if(sp.library_type == UNSTRANDED && ht.xs == '.') bb1->add_hit_intervals(ht, b1t);
		if(sp.library_type == UNSTRANDED && ht.xs == '.') bb2->add_hit_intervals(ht, b1t);
		*/
	}

	if(cfg.single_sample_multiple_threading == false) this->generate(bb1, mylock, index);
	else boost::asio::post(pool, [this, &mylock, bb1, index]{ this->generate(bb1, mylock, index); });
	index++;

	if(cfg.single_sample_multiple_threading == false) this->generate(bb2, mylock, index);
	else boost::asio::post(pool, [this, &mylock, bb2, index]{ this->generate(bb2, mylock, index); });
	index++;

	pool.join();

	return 0;
}

int generator::generate(bundle *bb, mutex &mylock, int index)
{
	if(bb == NULL) return 0;

	if(bb->tid < 0 || bb->hits.size() < cfg.min_num_hits_in_bundle) 
	{
		bb->clear();
		delete bb;
		return 0;
	}

	char buf[1024];
	strcpy(buf, hdr->target_name[bb->tid]);
	bb->chrm = string(buf);
	//bb->compute_strand(sp.library_type);

	//bb->print(index);

	splice_graph gr;
	graph_builder gb(*bb, cfg);
	gb.build(gr);
	gr.extend_strands();
	gr.build_vertex_index();

	//gr.print();
	//printf("\n");

	revise_splice_graph_full(gr, cfg);
	if(gr.count_junctions() <= 0) 
	{
		bb->clear();
		delete bb;
		return 0;
	}

	vector<pereads_cluster> vc;
	phase_set ps;

	graph_cluster gc(gr, bb->hits, cfg.max_reads_partition_gap, false);
	gc.build_pereads_clusters(vc);
	gc.build_phase_set_from_unpaired_reads(ps);

	bridge_solver bs(gr, vc, cfg, sp.insertsize_low, sp.insertsize_high);
	bs.build_phase_set(ps);
	//bs.print();

	vector<pereads_cluster> ub;
	bs.collect_unbridged_clusters(ub); 

	vector<splice_graph> grv;
	vector<phase_set> hsv;
	vector< vector<pereads_cluster> > ubv;
	partition(gr, ps, ub, grv, hsv, ubv);

	assert(grv.size() == hsv.size());
	assert(grv.size() == ubv.size());

	vector<combined_graph> tmp;
	for(int k = 0; k < grv.size(); k++)
	{
		if(grv[k].count_junctions() <= 0) continue;

		// TODO
		/*
		if(grv[k].mixed_strand_graph() == false) continue;
		grv[k].print();
		printf("\n");
		*/

		string gid = "gene." + tostring(index) + "." + tostring(k);
		combined_graph cb;
		cb.sid = sp.sample_id;
		cb.gid = gid;
		cb.build(grv[k], hsv[k], ubv[k]);

		//cb.print(k);
		//printf("\n");

		tmp.push_back(std::move(cb));
	}

	mylock.lock();
	for(int k = 0; k < tmp.size(); k++)
	{
		vcb.push_back(std::move(tmp[k]));
	}
	mylock.unlock();

	bb->clear();
	delete bb;
	return 0;
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
