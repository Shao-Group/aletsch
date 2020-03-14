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
	hit_base::operator=(h);
	itvm = h.itvm;
	itvi = h.itvi;
	itvd = h.itvd;
	return *this;
}

hit::hit(const hit &h)
	:hit_base(h)
{
	itvm = h.itvm;
	itvi = h.itvi;
	itvd = h.itvd;
}

hit::~hit()
{
}

hit::hit(bam1_t *b)
	:hit_base(b)
{} 

int hit::set_intervals(bam1_t *b)
{
	itvm.clear();
	itvi.clear();
	itvd.clear();
	int32_t p = pos;

	uint32_t *cigar = bam_get_cigar(b);

    for(int k = 0; k < n_cigar; k++)
	{
		if (bam_cigar_type(bam_cigar_op(cigar[k]))&2)
		{
			p += bam_cigar_oplen(cigar[k]);
		}

		if(bam_cigar_op(cigar[k]) == BAM_CMATCH)
		{
			int32_t s = p - bam_cigar_oplen(cigar[k]);
			itvm.push_back(pack(s, p));
		}

		if(bam_cigar_op(cigar[k]) == BAM_CINS)
		{
			itvi.push_back(pack(p - 1, p + 1));
		}

		if(bam_cigar_op(cigar[k]) == BAM_CDEL)
		{
			int32_t s = p - bam_cigar_oplen(cigar[k]);
			itvd.push_back(pack(s, p));
		}
	}
	return 0;
}
