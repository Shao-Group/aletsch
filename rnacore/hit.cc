/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#include <cstring>
#include <cassert>
#include <cstdio>
#include <sstream>
#include <cmath>

#include "hit.h"
#include "util.h"
#include "constants.h"

hit& hit::operator=(const hit &h)
{
	bam1_core_t::operator=(h);
	hid = h.hid;
	rpos = h.rpos;
	qname = h.qname;
	strand = h.strand;
	//spos = h.spos;
	xs = h.xs;
	ts = h.ts;
	hi = h.hi;
	nm = h.nm;
	nh = h.nh;
	return *this;
}

hit::hit(const hit &h)
	:bam1_core_t(h)
{
	hid = h.hid;
	rpos = h.rpos;
	qname = h.qname;
	strand = h.strand;
	//spos = h.spos;
	xs = h.xs;
	ts = h.ts;
	hi = h.hi;
	nm = h.nm;
	nh = h.nh;
}

hit::~hit()
{
}

hit::hit(bam1_t *b, int id)
	:bam1_core_t(b->core), hid(id)
{
	// fetch query name
	char buf[1024];
	char *qs = bam_get_qname(b);
	int l = strlen(qs);
	memcpy(buf, qs, l);
	buf[l] = '\0';
	qname = string(buf);

	// compute rpos
	rpos = pos + (int32_t)bam_cigar2rlen(n_cigar, bam_get_cigar(b));
}

bool hit::contain_splices(bam1_t *b) const
{
	uint32_t *cigar = bam_get_cigar(b);
    for(int k = 0; k < n_cigar; k++)
	{
		if(bam_cigar_op(cigar[k]) == BAM_CREF_SKIP) return true;
	}
	return false;
}

vector<int32_t> hit::extract_splices(bam1_t *b) const
{
	vector<int32_t> spos;
	uint32_t *cigar = bam_get_cigar(b);
	int32_t p = pos;
	int32_t q = 0;
    for(int k = 0; k < n_cigar; k++)
	{
		if (bam_cigar_type(bam_cigar_op(cigar[k]))&2)
			p += bam_cigar_oplen(cigar[k]);

		if (bam_cigar_type(bam_cigar_op(cigar[k]))&1)
			q += bam_cigar_oplen(cigar[k]);

		if(k == 0 || k == n_cigar - 1) continue;
		if(bam_cigar_op(cigar[k]) != BAM_CREF_SKIP) continue;
		//if(bam_cigar_op(cigar[k-1]) != BAM_CMATCH) continue;
		//if(bam_cigar_op(cigar[k+1]) != BAM_CMATCH) continue;
		////if(bam_cigar_oplen(cigar[k-1]) < min_flank) continue;
		////if(bam_cigar_oplen(cigar[k+1]) < min_flank) continue;

		int32_t s = p - bam_cigar_oplen(cigar[k]);
		//spos.push_back(pack(s, p));
		spos.push_back(s);
		spos.push_back(p);
	}
	return spos;
}

int hit::set_tags(bam1_t *b)
{
	ts = '.';
	uint8_t *p0 = bam_aux_get(b, "ts");
	if(p0 && (*p0) == 'A') ts = bam_aux2A(p0);

	xs = '.';
	uint8_t *p1 = bam_aux_get(b, "XS");
	if(p1 && (*p1) == 'A') xs = bam_aux2A(p1);

	if(xs == '.' && ts != '.')
	{
		// convert ts to xs
		if((flag & 0x10) >= 1 && ts == '+') xs = '-';
		if((flag & 0x10) >= 1 && ts == '-') xs = '+';
		if((flag & 0x10) <= 0 && ts == '+') xs = '+';
		if((flag & 0x10) <= 0 && ts == '-') xs = '-';
	}

	hi = -1;
	uint8_t *p2 = bam_aux_get(b, "HI");
	if(p2) hi = bam_aux2i(p2);

	nh = -1;
	uint8_t *p3 = bam_aux_get(b, "NH");
	if(p3) nh = bam_aux2i(p3);

	nm = 0;
	uint8_t *p4 = bam_aux_get(b, "nM");
	if(p4) nm = bam_aux2i(p4);

	uint8_t *p5 = bam_aux_get(b, "NM");
	if(p5) nm = bam_aux2i(p5);

	return 0;
}

bool hit::get_concordance() const
{
	if((flag & 0x10) <= 0 && (flag & 0x20) >= 1 && (flag & 0x40) >= 1 && (flag & 0x80) <= 0) return true;		// F1R2
	if((flag & 0x10) >= 1 && (flag & 0x20) <= 0 && (flag & 0x40) >= 1 && (flag & 0x80) <= 0) return true;		// R1F2
	if((flag & 0x10) <= 0 && (flag & 0x20) >= 1 && (flag & 0x40) <= 0 && (flag & 0x80) >= 1) return true;		// F2R1
	if((flag & 0x10) >= 1 && (flag & 0x20) <= 0 && (flag & 0x40) <= 0 && (flag & 0x80) >= 1) return true;		// R2F1
	return false;
}

int hit::set_strand(int libtype)
{
	strand = '.';
	
	if(libtype == FR_FIRST && ((flag & 0x1) >= 1))
	{
		if((flag & 0x10) <= 0 && (flag & 0x40) >= 1 && (flag & 0x80) <= 0) strand = '-';
		if((flag & 0x10) >= 1 && (flag & 0x40) >= 1 && (flag & 0x80) <= 0) strand = '+';
		if((flag & 0x10) <= 0 && (flag & 0x40) <= 0 && (flag & 0x80) >= 1) strand = '+';
		if((flag & 0x10) >= 1 && (flag & 0x40) <= 0 && (flag & 0x80) >= 1) strand = '-';
	}

	if(libtype == FR_SECOND && ((flag & 0x1) >= 1))
	{
		if((flag & 0x10) <= 0 && (flag & 0x40) >= 1 && (flag & 0x80) <= 0) strand = '+';
		if((flag & 0x10) >= 1 && (flag & 0x40) >= 1 && (flag & 0x80) <= 0) strand = '-';
		if((flag & 0x10) <= 0 && (flag & 0x40) <= 0 && (flag & 0x80) >= 1) strand = '-';
		if((flag & 0x10) >= 1 && (flag & 0x40) <= 0 && (flag & 0x80) >= 1) strand = '+';
	}

	if(libtype == FR_FIRST && ((flag & 0x1) <= 0))
	{
		if((flag & 0x10) <= 0) strand = '-';
		if((flag & 0x10) >= 1) strand = '+';
	}

	if(libtype == FR_SECOND && ((flag & 0x1) <= 0))
	{
		if((flag & 0x10) <= 0) strand = '+';
		if((flag & 0x10) >= 1) strand = '-';
	}

	return 0;
}

bool hit::operator<(const hit &h) const
{
	if(qname < h.qname) return true;
	if(qname > h.qname) return false;
	if(hi != -1 && h.hi != -1 && hi < h.hi) return true;
	if(hi != -1 && h.hi != -1 && hi > h.hi) return false;
	return (pos < h.pos);
}

int hit::print() const
{
	// print basic information
	printf("Hit %s: tid = %d, hid = %d, [%d-%d), mpos = %d, flag = %d, quality = %d, strand = %c, xs = %c, ts = %c, isize = %d, hi = %d\n", 
			qname.c_str(), tid, hid, pos, rpos, mpos, flag, qual, strand, xs, ts, isize, hi);

	return 0;

	/*
	printf(" start position (%d - )\n", pos);
	for(int i = 0; i < spos.size() / 2; i++)
	{
		int32_t p1 = spos[i * 2 + 0];
		int32_t p2 = spos[i * 2 + 1];
		printf(" splice position (%d - %d)\n", p1, p2);
	}
	printf(" end position (%d - )\n", rpos);
	return 0;
	*/
}

size_t hit::get_qhash() const
{
	return string_hash(qname);
}

/*
int hit::get_aligned_intervals(vector<int64_t> &v) const
{
	v.clear();
	int32_t p1 = pos;
	for(int k = 0; k < spos.size() / 2; k++)
	{
		int32_t p2 = spos[k * 2 + 0];
		v.push_back(pack(p1, p2));
		p1 = spos[k * 2 + 1];
	}
	v.push_back(pack(p1, rpos));
	return 0;
}

size_t hit::get_phash() const
{
	vector<int32_t> v;
	v.push_back(pos);
	v.insert(v.end(), spos.begin(), spos.end());
	v.push_back(rpos);
	return vector_hash(v);
}

bool hit::equal(const hit &h) const
{
	if(strand != h.strand) return false;

	if(tid != h.tid) return false;
	if(mtid != h.mtid) return false;
	if(pos != h.pos) return false;
	if(mpos != h.mpos) return false;
	if(rpos != h.rpos) return false;
	if(spos != h.spos) return false;
	if(flag != h.flag) return false;
	if(isize != h.isize) return false;
	if(n_cigar != h.n_cigar) return false;
	if(l_qname != h.l_qname) return false;
	if(l_extranul != h.l_extranul) return false;
	return true;
}
*/
