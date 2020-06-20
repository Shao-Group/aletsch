/*
Part of aletsch
(c) 2020 by Mingfu Shao, The Pennsylvania State University
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#include <cstdio>
#include <cassert>
#include <sstream>

#include "bundle.h"
#include "previewer.h"
#include "constants.h"
#include "graph_cluster.h"
#include "graph_builder.h"
#include "essential.h"

previewer::previewer(const parameters &c, sample_profile &s)
	: cfg(c), sp(s)
{
}

previewer::~previewer()
{
}

int previewer::infer_library_type()
{
	int total = 0;
	int single = 0;
	int paired = 0;

	int first = 0;
	int second = 0;
	vector<int> spn1;
	vector<int> spn2;

	int hid = 0;
	bam1_t *b1t = bam_init1();
	sp.open_align_file();

    while(sam_read1(sp.sfn, sp.hdr, b1t) >= 0)
	{
		if(total >= cfg.max_preview_reads) break;
		if(spn1.size() >= cfg.max_preview_spliced_reads && spn2.size() >= cfg.max_preview_spliced_reads) break;

		bam1_core_t &p = b1t->core;

		if((p.flag & 0x4) >= 1) continue;										// read is not mapped
		if((p.flag & 0x100) >= 1) continue;	// secondary alignment
		if(p.n_cigar > cfg.max_num_cigar) continue;									// ignore hits with more than max-num-cigar types
		if(p.qual < cfg.min_mapping_quality) continue;								// ignore hits with small quality
		if(p.n_cigar < 1) continue;												// should never happen

		total++;

		hit ht(b1t, hid++);
		ht.set_splices(b1t);
		ht.set_tags(b1t);

		if((ht.flag & 0x1) >= 1) paired ++;
		if((ht.flag & 0x1) <= 0) single ++;

		if(ht.xs == '.') continue;
		if(ht.xs == '+' && spn1.size() >= cfg.max_preview_spliced_reads) continue;
		if(ht.xs == '-' && spn2.size() >= cfg.max_preview_spliced_reads) continue;

		// predicted strand
		char xs = '.';

		// for paired read
		if((ht.flag & 0x1) >= 1 && (ht.flag & 0x10) <= 0 && (ht.flag & 0x20) >= 1 && (ht.flag & 0x40) >= 1 && (ht.flag & 0x80) <= 0) xs = '-';
		if((ht.flag & 0x1) >= 1 && (ht.flag & 0x10) >= 1 && (ht.flag & 0x20) <= 0 && (ht.flag & 0x40) <= 0 && (ht.flag & 0x80) >= 1) xs = '-';
		if((ht.flag & 0x1) >= 1 && (ht.flag & 0x10) >= 1 && (ht.flag & 0x20) <= 0 && (ht.flag & 0x40) >= 1 && (ht.flag & 0x80) <= 0) xs = '+';
		if((ht.flag & 0x1) >= 1 && (ht.flag & 0x10) <= 0 && (ht.flag & 0x20) >= 1 && (ht.flag & 0x40) <= 0 && (ht.flag & 0x80) >= 1) xs = '+';

		// for single read
		if((ht.flag & 0x1) <= 0 && (ht.flag & 0x10) <= 0) xs = '-';
		if((ht.flag & 0x1) <= 0 && (ht.flag & 0x10) >= 1) xs = '+';

		if(xs == '+' && xs == ht.xs) spn1.push_back(1);
		if(xs == '-' && xs == ht.xs) spn2.push_back(1);
		if(xs == '+' && xs != ht.xs) spn1.push_back(2);
		if(xs == '-' && xs != ht.xs) spn2.push_back(2);
	}

    bam_destroy1(b1t);
	sp.close_align_file();

	int spn = spn1.size() < spn2.size() ? spn1.size() : spn2.size();

	for(int k = 0; k < spn; k++)
	{
		if(spn1[k] == 1) first++;
		if(spn2[k] == 1) first++;
		if(spn1[k] == 2) second++;
		if(spn2[k] == 2) second++;
	}

	vector<string> vv;
	vv.push_back("empty");
	vv.push_back("unstranded");
	vv.push_back("first");
	vv.push_back("second");

	int s1 = UNSTRANDED;
	if(spn >= cfg.min_preview_spliced_reads && first > cfg.preview_infer_ratio * 2.0 * spn) s1 = FR_FIRST;
	if(spn >= cfg.min_preview_spliced_reads && second > cfg.preview_infer_ratio * 2.0 * spn) s1 = FR_SECOND;

	printf("infer-library-type (%s): reads = %d, single = %d, paired = %d, spliced = %d, first = %d, second = %d, inferred = %s\n",
			sp.align_file.c_str(), total, single, paired, spn, first, second, vv[s1 + 1].c_str());

	sp.library_type = s1;
	return 0;
}

int previewer::infer_insertsize()
{
	bundle bb1;
	bundle bb2;
	bb1.strand = '+';
	bb2.strand = '-';
	map<int32_t, int> m;
	int cnt = 0;
	int hid = 0;

	sp.open_align_file();
	bam1_t *b1t = bam_init1();
    while(sam_read1(sp.sfn, sp.hdr, b1t) >= 0)
	{
		bam1_core_t &p = b1t->core;

		if((p.flag & 0x4) >= 1) continue;										// read is not mapped
		if((p.flag & 0x100) >= 1) continue;	// secondary alignment
		if(p.n_cigar > cfg.max_num_cigar) continue;									// ignore hits with more than max-num-cigar types
		if(p.qual < cfg.min_mapping_quality) continue;								// ignore hits with small quality
		if(p.n_cigar < 1) continue;												// should never happen

		hit ht(b1t, hid++);
		ht.set_splices(b1t);
		ht.set_tags(b1t);
		ht.set_strand(sp.library_type);

		// truncate
		if(ht.tid != bb1.tid || ht.pos > bb1.rpos + cfg.min_bundle_gap)
		{
			cnt += process(bb1, m);
			bb1.clear();
			bb1.strand = '+';
		}
		if(ht.tid != bb2.tid || ht.pos > bb2.rpos + cfg.min_bundle_gap)
		{
			cnt += process(bb2, m);
			bb2.clear();
			bb2.strand = '-';
		}

		if(cnt >= 200000) break;

		// add hit
		if(cfg.uniquely_mapped_only == true && ht.nh != 1) continue;
		if(sp.library_type != UNSTRANDED && ht.strand == '+' && ht.xs == '-') continue;
		if(sp.library_type != UNSTRANDED && ht.strand == '-' && ht.xs == '+') continue;
		if(sp.library_type != UNSTRANDED && ht.strand == '.' && ht.xs != '.') ht.strand = ht.xs;
		if(sp.library_type != UNSTRANDED && ht.strand == '+') bb1.add_hit_intervals(ht, b1t);
		if(sp.library_type != UNSTRANDED && ht.strand == '-') bb2.add_hit_intervals(ht, b1t);
		if(sp.library_type == UNSTRANDED && ht.xs == '.') bb1.add_hit_intervals(ht, b1t);
		if(sp.library_type == UNSTRANDED && ht.xs == '.') bb2.add_hit_intervals(ht, b1t);
		if(sp.library_type == UNSTRANDED && ht.xs == '+') bb1.add_hit_intervals(ht, b1t);
		if(sp.library_type == UNSTRANDED && ht.xs == '-') bb2.add_hit_intervals(ht, b1t);
	}

    bam_destroy1(b1t);
	sp.close_align_file();

	int total = 0;
	for(map<int, int>::iterator it = m.begin(); it != m.end(); it++)
	{
		total += it->second;
	}

	if(total < 10000)
	{
		printf("not enough paired-end reads to create the profile (%d collected)\n", total);
		return 0;
	}

	int n = 0;
	double sx2 = 0;
	sp.insertsize_ave = 0;
	sp.insertsize_low = -1;
	sp.insertsize_high = -1;
	sp.insertsize_median = -1;
	vector<PI> vv(m.begin(), m.end());
	sort(vv.begin(), vv.end());
	//for(map<int32_t, int>::iterator it = m.begin(); it != m.end(); it++)
	for(int k = 0; k < vv.size(); k++)
	{
		n += vv[k].second;
		if(n >= 0.5 * total && sp.insertsize_median < 0) sp.insertsize_median = vv[k].first;
		sp.insertsize_ave += vv[k].second * vv[k].first;
		sx2 += vv[k].second * vv[k].first * vv[k].first;
		if(sp.insertsize_low == -1 && n >= 0.005 * total) sp.insertsize_low = vv[k].first;
		if(sp.insertsize_high == -1 && n >= 0.99 * total) sp.insertsize_high = vv[k].first;
		if(n >= 0.998 * total) break;
	}
	
	sp.insertsize_ave = sp.insertsize_ave * 1.0 / n;
	sp.insertsize_std = sqrt((sx2 - n * sp.insertsize_ave * sp.insertsize_ave) * 1.0 / n);

	printf("preview (%s) insertsize: sampled reads = %d, isize = %.2lf +/- %.2lf, median = %d, low = %d, high = %d\n", 
				sp.align_file.c_str(), total, sp.insertsize_ave, sp.insertsize_std, sp.insertsize_median, sp.insertsize_low, sp.insertsize_high);

	return 0;
}

int previewer::process(bundle &bd, map<int32_t, int> &m)
{
	if(bd.hits.size() < cfg.min_num_hits_in_bundle) return 0;
	if(bd.hits.size() > 20000) return 0;
	if(bd.tid < 0) return 0;

	char buf[1024];
	strcpy(buf, sp.hdr->target_name[bd.tid]);
	bd.chrm = string(buf);

	splice_graph gr;
	graph_builder gb(bd, cfg);
	gb.build(gr);
	//gr.print();

	gr.build_vertex_index();

	vector<pereads_cluster> vc;
	graph_cluster gc(gr, bd.hits, 2, false);
	gc.build_pereads_clusters(vc);
	
	//printf("process bundle in previewer with %lu hits, %lu clusters\n", bd.hits.size(), vc.size());

	int cnt = 0;
	for(int k = 0; k < vc.size(); k++)
	{
		pereads_cluster &pc = vc[k];
		int32_t p1 = pc.extend[1];
		int32_t p2 = pc.extend[2];
		int k1 = gr.locate_rbound(p1);
		int k2 = gr.locate_lbound(p2);

		if(k1 < 0 || k2 < 0 || k1 < k2) continue;

		vector<int32_t> chain;
		bool b = merge_intron_chains(pc.chain1, pc.chain2, chain);
		assert(b == true);

		int32_t length = get_total_length_of_introns(chain);
		int32_t d = pc.bounds[3] - pc.bounds[0] - length;

		//pc.print(k);
		//printf("d = %d, length = %d\n", d, length);

		cnt++;

		if(m.find(d) != m.end()) m[d]++;
		else m.insert(pair<int32_t, int>(d, 1));
		if(cnt >= 1000) return cnt;
	}
	return cnt;
}
