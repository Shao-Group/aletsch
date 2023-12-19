/*
Part of aletsch
(c) 2020 by Mingfu Shao, The Pennsylvania State University
See LICENSE for licensing.
*/

#include <cassert>
#include "transcript_set.h"
#include "constants.h"
//#include <boost/asio/post.hpp>
//#include <boost/asio/thread_pool.hpp>

trans_item::trans_item()
{}

trans_item::trans_item(const transcript &t, int c, int s)
{
	trst = t;
	count = c;
    trst.meta_tid = trst.transcript_id;
	if(samples.find(s) == samples.end()) 
    {
        samples.insert(make_pair(s,t));
        for(auto &x : samples) x.second.count2 = samples.size();
    }
    else 
    {
        samples[s].coverage = max(samples[s].coverage, t.coverage);
		//samples[s].alt_cov2 = t.cov2;
        //samples[s].cov2 = max(samples[s].cov2, t.cov2);
        samples[s].conf = max(samples[s].conf, t.conf);
        samples[s].abd = max(samples[s].abd, t.abd);
        samples[s].count1 = max(samples[s].count1, t.count1);
        printf("Check whether running!!!\n");
    }
}

int trans_item::merge(const trans_item &ti, int mode)
{
	if(mode == TRANSCRIPT_COUNT_ADD_COVERAGE_ADD) 
	{
		if(trst.exons.size() >= 2) trst.coverage += ti.trst.coverage;
		else if(trst.coverage < ti.trst.coverage) trst.coverage = ti.trst.coverage;

        trst.extend_bounds(ti.trst);
		count += ti.count;

        trst.cov2 = max(trst.cov2, ti.trst.cov2);
        //trst.alt_cov2 = ti.trst.cov2;
        trst.conf = max(trst.conf, ti.trst.conf);
        trst.abd = max(trst.abd, ti.trst.abd);
        trst.count1 = max(trst.count1, ti.trst.count1);

		for(auto &x : ti.samples)
		{
			if(samples.find(x.first) == samples.end()) samples.insert(x);
			else 
            {
                samples[x.first].cov2 = max(samples[x.first].cov2, x.second.cov2);
                //samples[x.first].alt_cov2 = x.second.cov2;
                samples[x.first].conf = max(samples[x.first].conf, x.second.conf);
                samples[x.first].abd = max(samples[x.first].abd, x.second.abd);
                samples[x.first].count1 = max(samples[x.first].count1, x.second.count1);
            }
		}

        trst.count2 = samples.size();
        for(auto &x : samples) 
        {
            x.second.coverage = trst.coverage;
            x.second.count2 = samples.size();
            x.second.meta_tid = trst.transcript_id;
        }
	}
	else if(mode == TRANSCRIPT_COUNT_ADD_COVERAGE_NUL) 
	{
		count += ti.count;
	}
	else assert(false);
	return 0;
}

int merge_sorted_trans_items(vector<trans_item> &vx, const vector<trans_item> &vy, int mode, double single_exon_ratio)
{
	vector<trans_item> vz;
	int kx = 0, ky = 0;
	while(kx < vx.size() && ky < vy.size())
	{
		int b = vx[kx].trst.compare1(vy[ky].trst, single_exon_ratio);
		if(b == 0)
		{
			vx[kx].merge(vy[ky], mode);
			vz.push_back(vx[kx]);
			kx++;
			ky++;
		}
		else if(b == 1)
		{
			vz.push_back(vx[kx]);
			kx++;
		}
		else if(b == -1)
		{
			vz.push_back(vy[ky]);
			ky++;
		}
		else assert(false);
	}

	assert(kx == vx.size() || ky == vy.size());

	for(int i = kx; i < vx.size(); i++) vz.push_back(vx[i]);
	for(int i = ky; i < vy.size(); i++) vz.push_back(vy[i]);

	vx.clear();
	vx = vz;
	return 0;
}

transcript_set::transcript_set(const string &c, int r, double s)
{
	chrm = c;
	rid = r;
	single_exon_overlap = s;
}

transcript_set::transcript_set(const transcript &t, int r, int count, int sid, double overlap)
{
	chrm = t.seqname;
	rid = r;
	single_exon_overlap = overlap;

	size_t h = t.get_intron_chain_hashing();
	trans_item ti(t, count, sid);
	vector<trans_item> v;
	v.push_back(std::move(ti));

	mt.insert(make_pair(h, v));
}

int transcript_set::clear()
{
	mt.clear();
	map<size_t, vector<trans_item>>().swap(mt);
	return 0;
}

int transcript_set::add(const transcript &t, int count, int sid, int mode)
{
	transcript_set ts(t, this->rid, count, sid, this->single_exon_overlap);
	add(ts, mode);
	return 0;
}

int transcript_set::add(const transcript_set &ts, int mode)
{
	if(ts.chrm != this->chrm) return 0;
	if(ts.rid != this->rid && this->rid != -9) return 0;

	for(auto &x : ts.mt)
	{
		map<size_t, vector<trans_item>>::iterator z = mt.find(x.first);
		if(z == mt.end())
		{
			mt.insert(x);
		}
		else
		{
			merge_sorted_trans_items(z->second, x.second, mode, single_exon_overlap);
		}
	}
	return 0;
}

int transcript_set::filter(int min_count)
{
	for(auto &x: mt)
	{
		vector<trans_item> v;
		for(auto &z: x.second)
		{
			if(z.count < min_count) continue;
			v.push_back(z);
		}
		x.second = std::move(v);
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
	printf("transcript-set: chrm = %s, rid = %d, mt.size() = %lu\n", chrm.c_str(), rid, mt.size());
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
		const transcript &x = v[k].trst;
		if(x.strand != t.strand) continue;
		bool b = x.equal1(t, single_exon_overlap);
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

int transcript_set_pool::clear()
{
	count = 0;
	for(int k = 0; k < tsets.size(); k++) tsets[k].clear();
	tsets.clear();
	vector<transcript_set>().swap(tsets);
	return 0;
}
