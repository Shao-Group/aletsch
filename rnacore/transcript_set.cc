#include "transcript_set.h"

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
			if(mode == 1 && z.coverage < t.coverage) z.coverage = t.coverage;
			if(mode == 2) z.coverage += t.coverage;
			z.count += t.count;
			found = true;
			break;
		}
		if(found == false) v.push_back(t);
	}
	return 0;
}

int transcript_set::add(transcript &t, int count, int mode)
{
	t.count = count;
	add(t, mode);
	return 0;
}

int transcript_set::add(const vector<transcript> &v, int mode)
{
	for(int i = 0; i < v.size(); i++)
	{
		add(v[i], mode);
	}
	return 0;
}

int transcript_set::add(const transcript_set &ts, int mode)
{
	for(auto &x : ts.mt)
	{
		add(x.second);
	}
	return 0;
}

int transcript_set::add_duplicates(const transcript_set &ts, int mode)
{
	vector<transcript> v = ts.get_duplicate_transcripts();
	add(v, mode);
	return 0;
}

vector<transcript> transcript_set::get_duplicate_transcripts() const
{
	vector<transcript> v;
	for(auto &x : mt)
	{
		for(auto &z : x.second)
		{
			if(z.count >= 2) v.push_back(z);
		}
	}
	return v;
}

vector<transcript> transcript_set::get_transcripts() const
{
	vector<transcript> v;
	for(auto &x : mt)
	{
		v.insert(v.end(), x.second.begin(), x.second.end());
	}
	return v;
}
