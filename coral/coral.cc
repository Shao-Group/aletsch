/*
Part of Coral package
(c) 2019 by Mingfu Shao, The Pennsylvania State University
See LICENSE for licensing.
*/

#include <cstdio>
#include <cassert>
#include <sstream>

#include "coral.h"
#include "bundle2.h"
#include "config.h"
#include "reference.h"
#include "htslib/bgzf.h"

coral::coral()
	: ref(ref_file)
{
    sfn = sam_open(input_file.c_str(), "r");
    hdr = sam_hdr_read(sfn);
    b1t = bam_init1();
	hid = 0;
	index = 0;
	if(output_bam_file != "")
	{
		bam_out = bgzf_open(output_bam_file.c_str(), "w");
		int f = bam_hdr_write(bam_out, hdr);
		if(f < 0) printf("fail to write header to %s\n", output_bam_file.c_str());
		if(f < 0) exit(0);
	}
}

coral::~coral()
{
    bam_destroy1(b1t);
    bam_hdr_destroy(hdr);
    sam_close(sfn);
	if(output_bam_file != "") bgzf_close(bam_out);
}

int coral::resolve()
{
	bb1.clear();
	bb2.clear();
	bb1.strand = '+';
	bb2.strand = '-';
    while(sam_read1(sfn, hdr, b1t) >= 0)
	{
		bam1_core_t &p = b1t->core;

		if(p.tid < 0) continue;
		if((p.flag & 0x4) >= 1) continue;										// read is not mapped
		if(p.n_cigar > max_num_cigar) continue;									// ignore hits with more than max-num-cigar types
		if(p.n_cigar < 1) continue;												// should never happen

		//if((p.flag & 0x100) >= 1 && use_second_alignment == false) continue;	// secondary alignment
		//if(p.qual < min_mapping_quality) continue;							// ignore hits with small quality

		hit ht(b1t, hid++);
		ht.set_splices(b1t, -1);
		ht.set_tags(b1t);
		ht.set_strand(library_type);

		//printf("tid = %d, "); ht.print();

		//if(ht.nh >= 2 && p.qual < min_mapping_quality) continue;
		//if(ht.nm > max_edit_distance) continue;

		// truncate
		bool b1 = false;
		if(ht.tid != bb1.tid || ht.pos > bb1.rpos + min_bundle_gap)
		{
			pool.push_back(bb1);
			bb1.clear();
			bb1.strand = '+';
			b1 = true;
		}

		bool b2 = false;
		if(ht.tid != bb2.tid || ht.pos > bb2.rpos + min_bundle_gap)
		{
			pool.push_back(bb2);
			bb2.clear();
			bb2.strand = '-';
			b2 = true;
		}

		// process
		if(b1 && b2) process();

		//printf("read strand = %c, xs = %c, ts = %c\n", ht.strand, ht.xs, ht.ts);

		// add hit
		if(uniquely_mapped_only == true && ht.nh != 1) continue;
		if(library_type != UNSTRANDED && ht.strand == '+' && ht.xs == '-') continue;
		if(library_type != UNSTRANDED && ht.strand == '-' && ht.xs == '+') continue;
		if(library_type != UNSTRANDED && ht.strand == '.' && ht.xs != '.') ht.strand = ht.xs;
		if(library_type != UNSTRANDED && ht.strand == '+') bb1.add_hit(ht);
		if(library_type != UNSTRANDED && ht.strand == '-') bb2.add_hit(ht);
		if(library_type == UNSTRANDED && ht.xs == '.') bb1.add_hit(ht);
		if(library_type == UNSTRANDED && ht.xs == '.') bb2.add_hit(ht);
		if(library_type == UNSTRANDED && ht.xs == '+') bb1.add_hit(ht);
		if(library_type == UNSTRANDED && ht.xs == '-') bb2.add_hit(ht);
	}

	pool.push_back(bb1);
	pool.push_back(bb2);
	process();

	write_multiple();
	return 0;
}

int coral::process()
{
	for(int i = 0; i < pool.size(); i++)
	{
		bundle_base &bb = pool[i];

		//printf("bundle %d has %lu reads\n", i, bb.hits.size());
		
		// TODO; do not use this parameter
		//if(bb.hits.size() < min_num_hits_in_bundle) continue;
		if(bb.tid < 0) continue;

		char buf[1024];
		strcpy(buf, hdr->target_name[bb.tid]);

		bundle2 bd(bb);
		bd.chrm = string(buf);
		bd.ref_trsts = ref.get_overlapped_transcripts(bd.chrm, bd.strand, bd.lpos, bd.rpos);
		bd.build();
		bd.print(index);
		if(output_bam_file != "") 
		{
			if(bd.strand == '+') bd.write_bam(bam_out, bmap1, umap1);
			if(bd.strand == '-') bd.write_bam(bam_out, bmap2, umap2);
		}
		index++;
	}

	write_unstrand();
	pool.clear();
	return 0;
}

int coral::write_unstrand()
{
	int n1 = 0;
	int n2 = 0;
	int n3 = 0;
	set<int> s;
	for(map<PI, bam1_t>::iterator it = bmap1.begin(); it != bmap1.end(); it++)
	{
		PI id = it->first;
		bam1_t &b = it->second;
		char c = '+';

		if(bmap2.find(id) == bmap2.end())
		{
			bam_aux_append(&b, "XS", 'A', 1, (uint8_t*)(&c));
			bam_write1(bam_out, &b);
			s.insert(id.first);
			s.insert(id.second);
			n1++;
		}
		assert(b.data != NULL);
		delete b.data;
	}

	for(map<PI, bam1_t>::iterator it = bmap2.begin(); it != bmap2.end(); it++)
	{
		PI id = it->first;
		bam1_t &b = it->second;
		char c = '-';

		if(bmap1.find(id) == bmap1.end())
		{
			bam_aux_append(&b, "XS", 'A', 1, (uint8_t*)(&c));
			bam_write1(bam_out, &b);
			s.insert(id.first);
			s.insert(id.second);
			n2++;
		}
		assert(b.data != NULL);
		delete b.data;
	}

	int u = 0;
	for(map<int, bam1_t>::iterator it = umap1.begin(); it != umap1.end(); it++)
	{
		int id = it->first;
		bam1_t &b = it->second;

		if(s.find(id) == s.end())
		{
			bam_write1(bam_out, &b);
			s.insert(id);
			u++;
		}

		assert(b.data != NULL);
		delete b.data;
	}

	for(map<int, bam1_t>::iterator it = umap2.begin(); it != umap2.end(); it++)
	{
		int id = it->first;
		bam1_t &b = it->second;

		if(s.find(id) == s.end())
		{
			bam_write1(bam_out, &b);
			s.insert(id);
			u++;
		}

		assert(b.data != NULL);
		delete b.data;
	}

	bmap1.clear();
	bmap2.clear();
	umap1.clear();
	umap2.clear();

	//printf("unstranded: bridged = (%d, %d, %d), unbridged = %d\n", n1, n2, n3, u);
	return 0;
}

int coral::write_multiple()
{
	int n = 0;
	for(int i = 0; i < uset.size(); i++)
	{
		bam1_t &b = uset[i];
		string s = hit::get_qname(&b);
		if(bset.find(s) == bset.end()) 
		{
			bam_write1(bam_out, &b);
			n++;
		}
		assert(b.data != NULL);
		delete b.data;
	}
	printf("total %lu bridged multiple-aligned reads; total %lu unbridged multiple-aligned reads: %d written\n",
			bset.size(), uset.size(), n);

	return 0;
}
