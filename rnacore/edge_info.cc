/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#include "edge_info.h"

edge_info::edge_info()
	: stddev(1.0), length(0)
{
	type = 0;
	jid = -1;
	strand = 0;
	weight = 0;
	count = 0;
    confidence = 0;
    abd = 0;
    samples.clear();
    spAbd.clear();
}

edge_info::edge_info(int l)
	: length(l)
{
	type = 0;
	jid = -1;
	count = 0;
	weight = 0;
	strand = 0;
    confidence = 0;
    abd = 0;
    samples.clear();
    spAbd.clear();
}

edge_info::edge_info(const edge_info &ei)
{
	stddev = ei.stddev;
	length = ei.length;
	type = ei.type;
	jid = ei.jid;
	count = ei.count;
	weight = ei.weight;
	strand = ei.strand;
    confidence = ei.confidence;
    abd = ei.abd;
    samples = ei.samples;
    spAbd = ei.spAbd;
}
