/*
Part of meta-scallop
(c) 2020 by Mingfu Shao, The Pennsylvania State University
See LICENSE for licensing.
*/

#include <cassert>
#include "transcript_set.h"
#include "constants.h"

int transcript_set::add(const trans_item &ti, int mode)
{
	transcript t = ti.trst;

	if(t.exons.size() <= 1) return 0;

	if(mt.size() == 0) strand = t.strand;
	if(mt.size() == 0) chrm = t.seqname;

	assert(t.strand == strand);
	assert(t.seqname == chrm);

	size_t h = t.get_intron_chain_hashing();

	auto it = mt.find(h);
	if(it == mt.end())
	{
		vector<trans_item> v;
		v.push_back(ti);
		mt.insert(make_pair(h, v));
	}
	else
	{
		auto &v = it->second;
		bool found = false;
		for(int k = 0; k < v.size(); k++)
		{
			bool b = v[k].trst.intron_chain_match(t);
			if(b == false) continue;

			if(mode == TRANSCRIPT_COUNT_ADD_COVERAGE_ADD) 
			{
				v[k].trst.coverage += t.coverage;
				v[k].trst.extend_bounds(t);
				v[k].count += ti.count;
				v[k].samples.insert(ti.samples.begin(), ti.samples.end());
			}
			else if(mode == TRANSCRIPT_COUNT_ADD_COVERAGE_NUL) 
			{
				v[k].count += ti.count;
			}
			else assert(false);

			found = true;
			break;
		}
		if(found == false) v.push_back(ti);
	}
	return 0;
}

int transcript_set::add(const transcript &t, int count, int sid, int mode)
{
	trans_item ti;
	ti.trst = t;
	ti.count = count;
	ti.samples.insert(sid);
	add(ti, mode);
	return 0;
}

int transcript_set::add(const transcript_set &ts, int min_count, int mode)
{
	for(auto &x : ts.mt)
	{
		for(auto &z : x.second)
		{
			if(z.count >= min_count) add(z, mode);
		}
	}
	return 0;
}

int transcript_set::increase_count(int count)
{
	for(auto &x : mt)
	{
		for(auto &z : x.second)
		{
			z.count += count;
		}
	}
	return 0;
}

int transcript_set::print() const
{
	printf("transcript-set: chrm = %s, strand = %c, mt.size() = %lu\n", chrm.c_str(), strand, mt.size());
	return 0;
}

vector<transcript> transcript_set::get_transcripts(int min_count) const
{
	vector<transcript> v;
	for(auto &x : mt)
	{
		for(auto &z : x.second)
		{
			if(z.count < min_count) continue;
			v.push_back(z.trst);
		}
	}
	return v;
}

pair<bool, trans_item> transcript_set::query(const transcript &t) const
{
	pair<bool, trans_item> p;
	size_t h = t.get_intron_chain_hashing();
	auto it = mt.find(h);
	if(it == mt.end()) 
	{
		p.first = false;
		return p;
	}

	auto &v = it->second;
	for(int k = 0; k < v.size(); k++)
	{
		transcript x = v[k].trst;
		if(x.strand != t.strand) continue;
		bool b = x.intron_chain_match(t);
		if(b == true) 
		{
			p.first = true;
			p.second = v[k];
			return p;
		}
	}
	p.first = true;
	return p;
}
