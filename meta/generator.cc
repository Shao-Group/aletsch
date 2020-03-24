/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#include <cstdio>
#include <cassert>
#include <sstream>
#include "boost/pending/disjoint_sets.hpp"

#include "constants.h"
#include "scallop/config.h"
#include "gtf.h"
#include "genome.h"
#include "generator.h"
#include "scallop.h"
#include "previewer.h"
#include "bridger.h"
#include "fcluster.h"
#include "essential.h"

generator::generator(vector<combined_graph> &v, const config &c)
	: vcb(v), cfg(c)
{
	previewer pre(cfg.input_file);
	pre.infer_library_type(cfg, sp);
	cfg.library_type = sp.library_type;
	pre.infer_insertsize(cfg, sp);

    sfn = sam_open(cfg.input_file.c_str(), "r");
    hdr = sam_hdr_read(sfn);
    b1t = bam_init1();
	index = 0;
	qlen = 0;
	qcnt = 0;
	//grout.open(gfile.c_str());
}

generator::~generator()
{
    bam_destroy1(b1t);
    bam_hdr_destroy(hdr);
    sam_close(sfn);
	//grout.close();
}

int generator::resolve()
{
    while(sam_read1(sfn, hdr, b1t) >= 0)
	{
		bam1_core_t &p = b1t->core;

		if((p.flag & 0x4) >= 1) continue;										// read is not mapped
		if((p.flag & 0x100) >= 1 && cfg.use_second_alignment == false) continue;	// secondary alignment
		if(p.n_cigar > cfg.max_num_cigar) continue;									// ignore hits with more than max-num-cigar types
		if(p.qual < cfg.min_mapping_quality) continue;								// ignore hits with small quality
		if(p.n_cigar < 1) continue;												// should never happen

		hit ht(b1t);
		ht.set_splices(b1t, cfg.min_flank_length);
		ht.set_tags(b1t);
		ht.set_strand(cfg.library_type);
		//ht.print();

		// TODO for test
		//if(ht.tid >= 1) break;

		//if(ht.nh >= 2 && p.qual < min_mapping_quality) continue;
		//if(ht.nm > max_edit_distance) continue;

		qlen += ht.qlen;
		qcnt += 1;

		// truncate
		if(ht.tid != bb1.tid || ht.pos > bb1.rpos + cfg.min_bundle_gap)
		{
			pool.push_back(bb1);
			bb1.clear();
		}
		if(ht.tid != bb2.tid || ht.pos > bb2.rpos + cfg.min_bundle_gap)
		{
			pool.push_back(bb2);
			bb2.clear();
		}

		// generate
		generate(cfg.batch_bundle_size);

		// add hit
		if(cfg.uniquely_mapped_only == true && ht.nh != 1) continue;
		if(cfg.library_type != UNSTRANDED && ht.strand == '+' && ht.xs == '-') continue;
		if(cfg.library_type != UNSTRANDED && ht.strand == '-' && ht.xs == '+') continue;
		if(cfg.library_type != UNSTRANDED && ht.strand == '.' && ht.xs != '.') ht.strand = ht.xs;
		if(cfg.library_type != UNSTRANDED && ht.strand == '+') bb1.add_hit_intervals(ht, b1t);
		if(cfg.library_type != UNSTRANDED && ht.strand == '-') bb2.add_hit_intervals(ht, b1t);
		if(cfg.library_type == UNSTRANDED && ht.xs == '.') bb1.add_hit_intervals(ht, b1t);
		if(cfg.library_type == UNSTRANDED && ht.xs == '.') bb2.add_hit_intervals(ht, b1t);
		if(cfg.library_type == UNSTRANDED && ht.xs == '+') bb1.add_hit_intervals(ht, b1t);
		if(cfg.library_type == UNSTRANDED && ht.xs == '-') bb2.add_hit_intervals(ht, b1t);
	}

	pool.push_back(bb1);
	pool.push_back(bb2);
	generate(0);

	return 0;
}

int generator::generate(int n)
{
	if(pool.size() < n) return 0;

	for(int i = 0; i < pool.size(); i++)
	{
		bundle_base &bb = pool[i];

		//printf("bundle %d has %lu reads\n", i, bb.hits.size());

		if(bb.hits.size() < cfg.min_num_hits_in_bundle) continue;
		if(bb.tid < 0) continue;

		char buf[1024];
		strcpy(buf, hdr->target_name[bb.tid]);

		bundle bd(bb, &cfg);

		bd.chrm = string(buf);
		bd.build();
		//bd.print(index);

		if(bd.gr.count_junctions() <= 0) continue;

		bridger br(bd.gr, bb.hits);
		br.length_median = sp.insertsize_median;
		br.length_low = sp.insertsize_low;
		br.length_high = sp.insertsize_high;
		br.resolve();

		//bd.gr.print();

		vector<fcluster> ub;
		br.collect_unbridged_fclusters(ub);

		//for(int k = 0; k < ub.size(); k++) ub[k].print(k);

		vector<splice_graph> grv;
		vector<hyper_set> hsv;
		vector< vector<fcluster> > ubv;
		partition(bd.gr, br.hs, ub, grv, hsv, ubv);
		
		assert(grv.size() == hsv.size());
		assert(grv.size() == ubv.size());

		for(int k = 0; k < grv.size(); k++)
		{
			if(grv[k].count_junctions() <= 0) continue;

			string gid = "gene." + tostring(index) + "." + tostring(k);
			combined_graph cb;
			cb.gid = gid;
			cb.build(grv[k], hsv[k], ubv[k]);

			//cb.print(k);	// TODO
			//printf("\n");

			vcb.push_back(cb);
		}

		//printf("\n");

		index++;
	}
	pool.clear();
	return 0;
}

int generator::partition(splice_graph &gr, hyper_set &hs, const vector<fcluster> &ub, vector<splice_graph> &grv, vector<hyper_set> &hsv, vector< vector<fcluster> > &ubv)
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

	// group with hyper_set
	for(MVII::const_iterator it = hs.nodes.begin(); it != hs.nodes.end(); it++)
	{
		const vector<int> &v = it->first;
		if(v.size() <= 1) continue;
		// here just verify all vertices in a phase are valid paths
		int p = ds.find_set(v[0]);
		for(int i = 1; i < v.size(); i++)
		{
			int q = ds.find_set(v[i]);
			if(p == q) continue;
			ds.link(p, q);
			p = ds.find_set(v[0]);
		}
	}

	// group with unbridged pairs
	for(int i = 0; i < ub.size(); i++)
	{
		if(ub[i].v1.size() <= 0) continue;
		if(ub[i].v2.size() <= 0) continue;
		int x = ub[i].v1.back();
		int y = ub[i].v2.front();
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

	hsv.resize(vv.size());
	for(MVII::const_iterator it = hs.nodes.begin(); it != hs.nodes.end(); it++)
	{
		vector<int> v = it->first;
		int c = it->second;
		if(v.size() <= 0) continue;
		int p = ds.find_set(v[0]);
		assert(m.find(p) != m.end());
		int k = m[p];
		assert(k >= 0 && k < vv.size());
		vector<int> vv = project_vector(v, vm[k]);
		assert(vv.size() == v.size());
		for(int i = 0; i < vv.size(); i++) vv[i]--;
		hsv[k].add_node_list(vv, c);
	}

	ubv.resize(vv.size());
	for(int i = 0; i < ub.size(); i++)
	{
		if(ub[i].v1.size() <= 0) continue;
		if(ub[i].v2.size() <= 0) continue;
		int x = ub[i].v1.back();
		int y = ub[i].v2.front();
		int p = ds.find_set(x);
		assert(p == ds.find_set(y));
		assert(m.find(p) != m.end());
		int k = m[p];
		assert(k >= 0 && k < vv.size());

		fcluster fc = ub[i];
		fc.v1 = project_vector(fc.v1, vm[k]);
		fc.v2 = project_vector(fc.v2, vm[k]);
		//assert(fc.v1.size() == ub[i].v1.size());
		//assert(fc.v2.size() == ub[i].v2.size());
		if(fc.v1.size() != ub[i].v1.size()) continue;
		if(fc.v2.size() != ub[i].v2.size()) continue;
		ubv[k].push_back(fc);
	}

	return 0;
}

int generator::write_graph(splice_graph &gr, hyper_set &hs)
{
	gr.write(grout);

	for(MVII::const_iterator it = hs.nodes.begin(); it != hs.nodes.end(); it++)
	{
		const vector<int> &v = it->first;
		int c = it->second;
		gr.write(grout, v, c, "phase");
	}

	scallop sc(gr, hs, &cfg);
	sc.preassemble();

	for(int i = 0; i < sc.paths.size(); i++)
	{
		const vector<int> &v = sc.paths[i].v;
		double w = sc.paths[i].abd;
		gr.write(grout, v, w, "path");
	}

	grout << endl;
	return 0;
}
