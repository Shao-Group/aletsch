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

int transcript_set::add(const transcript &t, int mode)
{
	if(t.exons.size() <= 1) return 0;
	assert(t.count >= 1);

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
			if(z.strand != t.strand) continue;
			bool b = z.intron_chain_match(t);
			if(b == false) continue;
			if(mode == ADD_TRANSCRIPT_COVERAGE_SUM) z.coverage += t.coverage;
			if(mode == ADD_TRANSCRIPT_COVERAGE_MAX && z.coverage < t.coverage) z.coverage = t.coverage;
			if(mode == ADD_TRANSCRIPT_COVERAGE_MIN && z.coverage > t.coverage) z.coverage = t.coverage;
			z.extend_bounds(t);
			found = true;
			break;
		}
		if(found == false) v.push_back(t);
	}
	return 0;
}

int transcript_set::add(const transcript_set &ts, int min_count, int mode)
{
	vector<transcript> v = ts.get_transcripts(min_count);
	for(int i = 0; i < v.size(); i++) add(v[i], mode);
	return 0;
}

int transcript_set::add(const transcript_set &ts, int mode)
{
	add(ts, 1, mode);
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

int transcript_set::print()
{
	vector<transcript> vv = get_transcripts(1);
	for(int i = 0; i < vv.size(); i++) vv[i].write(cout);
	return 0;
}
