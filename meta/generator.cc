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
#include "super_graph.h"
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

		// process
		process(cfg.batch_bundle_size);

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
	process(0);

	return 0;
}

int generator::process(int n)
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

		//generate(bd.gr, br.hs);

		vector<fcluster> ub;
		// TODO
		//br.collect_unbridged_fclusters(ub);

		vector< set<int> > vv;
		vector< vector<int> > uv;
		//partition(bd.gr, ub, vv, uv);
		assert(vv.size() == uv.size());

		//bd.gr.print();
		//br.hs.print_nodes();
		for(int k = 0; k < vv.size(); k++)
		{
			splice_graph gr;
			hyper_set hs;
			build_child_splice_graph(bd.gr, gr, vv[k]);
			build_child_hyper_set(br.hs, hs, vv[k]);

			if(gr.count_junctions() <= 0) continue;
			
			/*
			printf("---\n");
			printf("set = ( ");
			printv(vector<int>(vv[k].begin(), vv[k].end()));
			printf(")\n");
			gr.print();
			hs.print_nodes();
			*/

			vector<fcluster> vf;
			for(int j = 0; j < uv[k].size(); j++) vf.push_back(ub[uv[k][j]]);

			string gid = "gene." + tostring(index) + "." + tostring(k);
			combined_graph cb;
			cb.build(gr, hs, vf);
			vcb.push_back(cb);
		}

		//printf("\n");

		index++;
	}
	pool.clear();
	return 0;
}

int generator::partition(splice_graph &gr, const vector<fcluster> &ub, vector< set<int> > &vv, vector< vector<int> > &uv)
{
	int n = gr.num_vertices();

	vector<int> rank(n, -1);
	vector<int> parent(n, -1);

	disjoint_sets<int*, int*> ds(&rank[0], &parent[0]);
	for(int k = 0; k < n; k++) ds.make_set(k);

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

	for(int i = 0; i < ub.size(); i++)
	{
		if(ub[i].v1.size() <= 0) continue;
		if(ub[i].v2.size() <= 0) continue;
		int x = ub[i].v1.back();
		int y = ub[i].v2.front();
		ds.union_set(x, y);
	}

	map<int, int> m;
	vv.clear();
	uv.clear();
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

	uv.resize(vv.size());
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
		assert(k >= 0 && k < uv.size());
		uv[k].push_back(i);
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
