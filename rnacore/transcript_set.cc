#include <cassert>
#include "transcript_set.h"
#include "constants.h"

bool transcript_set::query(const transcript &t) const
{
	size_t h = t.get_intron_chain_hashing();
	auto it = mt.find(h);
	if(it == mt.end()) return false;

	auto &v = it->second;
	for(int k = 0; k < v.size(); k++)
	{
		if(v[k].strand != t.strand) continue;
		bool b = v[k].intron_chain_match(t);
		if(b == true) return true;
	}
	return false;
}

int transcript_set::add(const transcript &t1, int mode)
{
	transcript t = t1;

	if(t.exons.size() <= 1) return 0;
	assert(t.count >= 1);

	if(mt.size() == 0) strand = t.strand;
	if(mt.size() == 0) chrm = t.seqname;

	assert(t.strand == strand);
	assert(t.seqname == chrm);

	size_t h = t.get_intron_chain_hashing();

	auto it = mt.find(h);
	if(it == mt.end())
	{
		vector<transcript> v;
		v.push_back(t);
		mt.insert(make_pair(h, v));
	}
	else
	{
		auto &v = it->second;
		bool found = false;
		for(int k = 0; k < v.size(); k++)
		{
			transcript &z = v[k];
			bool b = z.intron_chain_match(t);
			if(b == false) continue;

			if(mode == TRANSCRIPT_COUNT_ADD_COVERAGE_ADD) 
			{
				z.count += t.count;
				z.coverage += t.coverage;
				z.extend_bounds(t);
			}

			if(mode == TRANSCRIPT_COUNT_ADD_COVERAGE_NUL) 
			{
				z.count += t.count;
			}

			if(mode == TRANSCRIPT_COUNT_MAX_COVERAGE_MAX)
			{
				if(t.count > z.count) z = t;
				if(t.count == z.count && t.coverage > z.coverage) z.coverage = t.coverage;
				if(t.count == z.count && t.coverage > z.coverage) z.extend_bounds(t);
			}

			found = true;
			break;
		}
		if(found == false) v.push_back(t);
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

int transcript_set::add(const transcript_set &ts, int mode)
{
	for(auto &x : ts.mt)
	{
		for(auto &z : x.second)
		{
			add(z, mode);
		}
	}
	return 0;
}

vector<transcript> transcript_set::get_transcripts(int min_count) const
{
	vector<transcript> v;
	for(auto &x : mt)
	{
		for(auto &z : x.second)
		{
			if(z.count >= min_count) v.push_back(z);
		}
	}
	return v;
}

int transcript_set::print() const
{
	printf("transcript-set: chrm = %s, strand = %c, mt.size() = %lu\n", chrm.c_str(), strand, mt.size());
	vector<transcript> vv = get_transcripts(1);
	for(int i = 0; i < vv.size(); i++) vv[i].write(cout);
	return 0;
}
