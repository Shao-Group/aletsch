#include "hit_core.h"

hit_core::hit_core(const hit_core &h)
	:bam1_core_t(h)
{
	rpos = h.rpos;
	qname = h.qname;
	hi = h.hi;
	nh = h.nh;
	xs = h.xs;
}

hit_core::hit_core(const hit &h)
	:bam1_core_t(h)
{
	rpos = h.rpos;
	qname = h.qname;
	hi = h.hi;
	nh = h.nh;
	xs = h.xs;
}
