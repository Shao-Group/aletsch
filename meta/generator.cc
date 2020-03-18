/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#include <cstdio>
#include <cassert>
#include <sstream>

#include "constants.h"
#include "scallop/config.h"
#include "gtf.h"
#include "genome.h"
#include "generator.h"
#include "scallop.h"
#include "super_graph.h"
#include "previewer.h"
#include "bridger.h"

generator::generator(vector<combined_graph> &v, const config &c)
	: vcb(v), cfg(c)
{
	previewer pre(cfg.input_file);
	pre.infer_library_type(cfg, sp);
	pre.infer_insertsize(cfg, sp);
	cfg.library_type = sp.library_type;

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

		bridger br(bd.gr, bb.hits);
		br.length_median = sp.insertsize_median;
		br.length_low = sp.insertsize_low;
		br.length_high = sp.insertsize_high;
		br.resolve();

		generate(bd.gr, br.hs);
		index++;
	}
	pool.clear();
	return 0;
}

int generator::generate(const splice_graph &gr0, const hyper_set &hs0)
{
	super_graph sg(gr0, hs0, &cfg);
	sg.build();

	vector<transcript> gv;
	for(int k = 0; k < sg.subs.size(); k++)
	{
		string gid = "gene." + tostring(index) + "." + tostring(k);

		if(cfg.verbose >= 2 && k == 0) sg.print();

		splice_graph &gr = sg.subs[k];
		hyper_set &hs = sg.hss[k];
		gr.gid = gid;

		if(gr.count_junctions() >= 1)
		{
			//write_graph(gr, hs);
			combined_graph cb;
			cb.build(gr, hs);
			vcb.push_back(cb);
		}
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
