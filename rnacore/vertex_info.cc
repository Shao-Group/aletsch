/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#include "vertex_info.h"

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
    supports.clear();
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
    supports.clear();
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
    supports = vi.supports;
}
