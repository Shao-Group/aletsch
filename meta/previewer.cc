/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#include <cstdio>
#include <cassert>
#include <sstream>

#include "previewer.h"
#include "constants.h"
#include "bridger.h"
#include "bundle.h"

previewer::previewer(const string &file)
{
	input_file = file;
	max_preview_reads = 2000000;
	max_preview_spliced_reads = 50000;
	min_preview_spliced_reads = 10000;
	preview_infer_ratio = 0.8;
}

previewer::~previewer()
{
}

int previewer::open_file()
{
    sfn = sam_open(input_file.c_str(), "r");
    hdr = sam_hdr_read(sfn);
    b1t = bam_init1();
	return 0;
}

int previewer::close_file()
{
    bam_destroy1(b1t);
    bam_hdr_destroy(hdr);
    sam_close(sfn);
	return 0;
}

int previewer::infer_library_type(config &cfg, sample_profile &sx)
{
	open_file();

	int total = 0;
	int single = 0;
	int paired = 0;

	int first = 0;
	int second = 0;
	vector<int> sp1;
	vector<int> sp2;

    while(sam_read1(sfn, hdr, b1t) >= 0)
	{
		if(total >= max_preview_reads) break;
		if(sp1.size() >= max_preview_spliced_reads && sp2.size() >= max_preview_spliced_reads) break;

		bam1_core_t &p = b1t->core;

		if((p.flag & 0x4) >= 1) continue;										// read is not mapped
		if((p.flag & 0x100) >= 1) continue;	// secondary alignment
		if(p.n_cigar > cfg.max_num_cigar) continue;									// ignore hits with more than max-num-cigar types
		if(p.qual < cfg.min_mapping_quality) continue;								// ignore hits with small quality
		if(p.n_cigar < 1) continue;												// should never happen

		total++;

		hit ht(b1t);
		ht.set_splices(b1t, cfg.min_flank_length);
		ht.set_tags(b1t);

		if((ht.flag & 0x1) >= 1) paired ++;
		if((ht.flag & 0x1) <= 0) single ++;

		if(ht.xs == '.') continue;
		if(ht.xs == '+' && sp1.size() >= max_preview_spliced_reads) continue;
		if(ht.xs == '-' && sp2.size() >= max_preview_spliced_reads) continue;

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

		if(xs == '+' && xs == ht.xs) sp1.push_back(1);
		if(xs == '-' && xs == ht.xs) sp2.push_back(1);
		if(xs == '+' && xs != ht.xs) sp1.push_back(2);
		if(xs == '-' && xs != ht.xs) sp2.push_back(2);
	}

	int sp = sp1.size() < sp2.size() ? sp1.size() : sp2.size();

	for(int k = 0; k < sp; k++)
	{
		if(sp1[k] == 1) first++;
		if(sp2[k] == 1) first++;
		if(sp1[k] == 2) second++;
		if(sp2[k] == 2) second++;
	}

	vector<string> vv;
	vv.push_back("empty");
	vv.push_back("unstranded");
	vv.push_back("first");
	vv.push_back("second");

	int s1 = UNSTRANDED;
	if(sp >= min_preview_spliced_reads && first > preview_infer_ratio * 2.0 * sp) s1 = FR_FIRST;
	if(sp >= min_preview_spliced_reads && second > preview_infer_ratio * 2.0 * sp) s1 = FR_SECOND;

	printf("infer-library-type (%s): reads = %d, single = %d, paired = %d, spliced = %d, first = %d, second = %d, inferred = %s\n",
			input_file.c_str(), total, single, paired, sp, first, second, vv[s1 + 1].c_str());

	sx.library_type = s1;
	close_file();
	return 0;
}

int previewer::infer_insertsize(config &cfg, sample_profile &sp)
{
	//printf("preview insertsize for file %s\n", input_file.c_str());
	open_file();

	bundle_base bb1;
	bundle_base bb2;
	bb1.strand = '+';
	bb2.strand = '-';
	map<int32_t, int> m;
	int cnt = 0;

    while(sam_read1(sfn, hdr, b1t) >= 0)
	{
		bam1_core_t &p = b1t->core;

		if((p.flag & 0x4) >= 1) continue;										// read is not mapped
		if((p.flag & 0x100) >= 1) continue;	// secondary alignment
		if(p.n_cigar > cfg.max_num_cigar) continue;									// ignore hits with more than max-num-cigar types
		if(p.qual < cfg.min_mapping_quality) continue;								// ignore hits with small quality
		if(p.n_cigar < 1) continue;												// should never happen

		hit ht(b1t);
		ht.set_splices(b1t, cfg.min_flank_length);
		ht.set_tags(b1t);
		ht.set_strand(cfg.library_type);

		// truncate
		if(ht.tid != bb1.tid || ht.pos > bb1.rpos + cfg.min_bundle_gap)
		{
			cnt += process(bb1, cfg, m);
			bb1.clear();
			bb1.strand = '+';
		}
		if(ht.tid != bb2.tid || ht.pos > bb2.rpos + cfg.min_bundle_gap)
		{
			cnt += process(bb2, cfg, m);
			bb2.clear();
			bb2.strand = '-';
		}

		if(cnt >= 200000) break;

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

	int total = 0;
	for(map<int, int>::iterator it = m.begin(); it != m.end(); it++)
	{
		total += it->second;
	}

	if(total < 10000)
	{
		printf("not enough paired-end reads to create the profile (%d collected)\n", total);
		close_file();
		return 0;
	}

	int n = 0;
	double sx2 = 0;
	sp.insertsize_ave = 0;
	sp.insertsize_low = -1;
	sp.insertsize_high = -1;
	sp.insertsize_median = -1;
	for(map<int32_t, int>::iterator it = m.begin(); it != m.end(); it++)
	{
		n += it->second;
		if(n >= 0.5 * total && sp.insertsize_median < 0) sp.insertsize_median = it->first;
		sp.insertsize_ave += it->second * it->first;
		sx2 += it->second * it->first * it->first;
		if(sp.insertsize_low == -1 && n >= 0.01 * total) sp.insertsize_low = it->first;
		if(sp.insertsize_high == -1 && n >= 0.98 * total) sp.insertsize_high = it->first;
		if(n >= 0.998 * total) break;
	}
	
	sp.insertsize_ave = sp.insertsize_ave * 1.0 / n;
	sp.insertsize_std = sqrt((sx2 - n * sp.insertsize_ave * sp.insertsize_ave) * 1.0 / n);

	/*
	sp.insertsize_profile.assign(sp.insertsize_high, 1);
	n = sp.insertsize_high;
	for(map<int32_t, int>::iterator it = m.begin(); it != m.end(); it++)
	{
		if(it->first >= sp.insertsize_high) continue;
		sp.insertsize_profile[it->first] += it->second;
		n += it->second;
	}

	for(int i = 0; i < sp.insertsize_profile.size(); i++)
	{
		sp.insertsize_profile[i] = sp.insertsize_profile[i] * 1.0 / n;
	}
	*/

	printf("preview (%s) insertsize: sampled reads = %d, isize = %.2lf +/- %.2lf, median = %d, low = %d, high = %d\n", 
				input_file.c_str(), total, sp.insertsize_ave, sp.insertsize_std, sp.insertsize_median, sp.insertsize_low, sp.insertsize_high);

	close_file();
	return 0;
}

int previewer::process(bundle_base &bb, config &cfg, map<int32_t, int> &m)
{
	if(bb.hits.size() < cfg.min_num_hits_in_bundle) return 0;
	if(bb.hits.size() > 20000) return 0;
	if(bb.tid < 0) return 0;

	char buf[1024];
	strcpy(buf, hdr->target_name[bb.tid]);

	bundle bd(bb, &cfg);

	bd.chrm = string(buf);
	bd.build();
	//bd.print(index);

	bridger br(bd.gr, bb.hits);
	br.length_low = 0;
	br.length_high = 9999;
	br.build_vertex_index();
	br.build_fragments();
	br.build_fclusters();
	br.vote();

	int cnt = 0;
	for(int k = 0; k < br.fclusters.size(); k++)
	{
		fcluster &fc = br.fclusters[k];
		if(fc.bbp.type != 1) continue;

		for(int j = 0; j < fc.frset.size(); j++)
		{
			fragment &fr = br.fragments[fc.frset[j]];
			int32_t len = br.compute_aligned_length(fr, fc.bbp.v);

			cnt++;

			if(m.find(len) != m.end()) m[len]++;
			else m.insert(pair<int32_t, int>(len, 1));
			if(cnt >= 1000) return cnt;
		}
	}
	return cnt;
}
