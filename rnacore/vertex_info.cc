/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#include "vertex_info.h"
#include <cfloat>

vertex_info::vertex_info()
{
	stddev = 1.0;
	maxcov = 0;
	length = 0;
	sdist = -1;
	tdist = -1;
	type = -1;
	lpos = 0;
	rpos = 0;
	pos = 0;
	count = 0;
	lstrand = '.';
	rstrand = '.';
	regional = false;

    trstSupport = 0;
	indel_sum_cov = 0;
    indel_ratio = 0;
    left_indel = -1;
    right_indel = -1;

    boundary_loss1 = 0;
    boundary_loss2 = 0;
    boundary_loss3 = 0;
    boundary_merged_loss = 0;
    unbridge_leaving_count = 0;
    unbridge_leaving_ratio = 0;
    unbridge_coming_count = 0;
    unbridge_coming_ratio = 0;
}

vertex_info::vertex_info(int l)
	: length(l)
{
	stddev = 1.0;
	maxcov = 0;
	sdist = -1;
	tdist = -1;
	type = -1;
	lpos = 0;
	rpos = 0;
	pos = 0;
	count = 0;
	lstrand = '.';
	rstrand = '.';
	regional = false;

	trstSupport = 0;
    indel_sum_cov = 0;
    indel_ratio = 0;
    left_indel = -1;
    right_indel = -1;

    boundary_loss1 = 0;
    boundary_loss2 = 0;
    boundary_loss3 = 0;
    boundary_merged_loss = 0;
    unbridge_leaving_count = 0;
    unbridge_leaving_ratio = 0;
    unbridge_coming_count = 0;
    unbridge_coming_ratio = 0;


}

vertex_info::vertex_info(const vertex_info &vi)
{
	stddev = vi.stddev;
	maxcov = vi.maxcov;
	length = vi.length;
	sdist = vi.sdist;
	tdist = vi.tdist;
	type = vi.type;
	lpos = vi.lpos;
	rpos = vi.rpos;
	pos = vi.pos;
	count = vi.count;
	lstrand = vi.lstrand;
	rstrand = vi.rstrand;
	regional = vi.regional;

	trstSupport = vi.trstSupport;
    indel_sum_cov = vi.indel_sum_cov;
    indel_ratio = vi.indel_ratio;
    left_indel = vi.left_indel;
    right_indel = vi.right_indel;

    boundary_loss1 = vi.boundary_loss1;
    boundary_loss2 = vi.boundary_loss2;
    boundary_loss3 = vi.boundary_loss3;
    boundary_merged_loss = vi.boundary_merged_loss;
    unbridge_leaving_count = vi.unbridge_leaving_count;
    unbridge_leaving_ratio = vi.unbridge_leaving_ratio;
    unbridge_coming_count = vi.unbridge_coming_count;
    unbridge_coming_ratio = vi.unbridge_coming_ratio;
}
