/*
Part of aletsch
(c) 2020 by Mingfu Shao, The Pennsylvania State University
See LICENSE for licensing.
*/

#include <algorithm>
#include <sstream>
#include <fstream>
#include "bundle_group.h"
#include "parameters.h"
#include <boost/asio/post.hpp>
#include <boost/asio/thread_pool.hpp>

//mutex bundle_group::gmutex;

bundle_group::bundle_group(string c, char s, int r, const parameters &f, thread_pool &p)
	: cfg(f), tpool(p), tmerge(c, r, f.min_single_exon_clustering_overlap)
{
	chrm = c;
	strand = s;
	rid = r;
}

int bundle_group::resolve()
{
	build_splices();
	build_splice_index();

	disjoint_set ds(gset.size());
	grouped.assign(gset.size(), false);

	// round one
	for(auto &z: sindex)
	{
		if(z.second.size() <= 1) continue;
		process_subset(z.second, ds, cfg.max_grouping_similarity);
	}

	stats(ds, 1);
	if(cfg.verbose >= 2) print();

	// round two
	for(auto &z: sindex)
	{
		if(z.second.size() <= 1) continue;
		process_subset(z.second, ds, cfg.min_grouping_similarity);
	}

	build_groups(ds);

	stats(ds, 2);
	if(cfg.verbose >= 2) print();

	return 0;
}

int bundle_group::clear()
{
	for(int k = 0; k < splices.size(); k++)
	{
		splices[k].clear();
		vector<int32_t>().swap(splices[k]);
	}
	splices.clear();
	vector<vector<int32_t>>().swap(splices);

	for(auto &z: jmaps) join_interval_map().swap(z);
	jmaps.clear();
	vector<join_interval_map>().swap(jmaps);

	for(int k = 0; k < gvv.size(); k++)
	{
		gvv[k].clear();
		vector<int>().swap(gvv[k]);
	}
	gvv.clear();
	vector<vector<int>>().swap(gvv);

	sindex.clear();
	MISI().swap(sindex);

	jindex.clear();
	interval_set_map().swap(jindex);

	grouped.clear();
	vector<bool>().swap(grouped);
	return 0;
}

int bundle_group::process_subset(const set<int> &s, disjoint_set &ds, double d)
{
	vector<int> ss;
	filter(s, ds, ss);

	vector<PPID> vpid;
	build_splice_similarity(ss, vpid, ds, false, d);

	augment_disjoint_set(vpid, ds);
	return 0;
}

int bundle_group::build_splices()
{
	splices.clear();
	for(int i = 0; i < gset.size(); i++)
	{
		vector<int32_t> v = gset[i].hcst.get_splices();
		splices.push_back(std::move(v));
	}
	return 0;
}

int bundle_group::build_join_interval_maps()
{
	jmaps.clear();
	for(int i = 0; i < gset.size(); i++)
	{
		join_interval_map jmap;
		for(SIMI it = gset[i].mmap.begin(); it != gset[i].mmap.end(); it++)
		{
			jmap += make_pair(it->first, 1);
		}
		jmaps.push_back(std::move(jmap));
	}
	return 0;
}

int bundle_group::build_splice_index()
{
	sindex.clear();
	for(int k = 0; k < gset.size(); k++)
	{
		for(int i = 0; i < splices[k].size(); i++)
		{
			int32_t p = splices[k][i];
			MISI::iterator it = sindex.find(p);
			if(it == sindex.end())
			{
				set<int> s;
				s.insert(k);
				sindex.insert(PISI(p, s));
			}
			else
			{
				it->second.insert(k);
			}
		}
	}
	return 0;
}

int bundle_group::build_join_interval_map_index()
{
	jindex.clear();
	for(int k = 0; k < jmaps.size(); k++)
	{
		if(jmaps[k].begin() == jmaps[k].end()) continue;
		int32_t x = lower(jmaps[k].begin()->first);
		int32_t y = upper(jmaps[k].rbegin()->first);

		set<int> s;
		s.insert(k);
		jindex += make_pair(interval32(x, y), s);
	}
	return 0;
}

int bundle_group::build_splice_similarity(const vector<int> &ss, vector<PPID> &vpid, disjoint_set &ds, bool local, double min_similarity)
{
	//printf("START BUILD SIMILARITY; cfg.max_num_junctions_to_combine = %d, ss.size() =%lu\n", cfg.max_num_junctions_to_combine, ss.size());
	for(int xi = 0; xi < ss.size(); xi++)
	{
		int i = ss[xi];
		if(splices[i].size() / 2.0 > cfg.max_num_junctions_to_combine) continue;
		int pi = ds.find_set(i);
		for(int xj = 0; xj < ss.size(); xj++)
		{
			int j = ss[xj];
			if(i >= j) continue;
			if(splices[j].size() / 2.0 > cfg.max_num_junctions_to_combine) continue;

			int pj = ds.find_set(j);
			if(pi == pj) continue;

			assert(gset[i].chrm == gset[j].chrm);
			assert(gset[i].strand == gset[j].strand);

			vector<int32_t> vv(splices[i].size() + splices[j].size(), 0);
			vector<int32_t>::iterator it = set_intersection(splices[i].begin(), splices[i].end(), splices[j].begin(), splices[j].end(), vv.begin());
			int c = it - vv.begin();
			//double r = c * 1.0 / (splices[i].size() + splices[j].size() - c);
			int small = splices[i].size() < splices[j].size() ? splices[i].size() : splices[j].size();
			double r = c * 1.0 / small;

			if(cfg.verbose >= 2) printf("graph-similarity: r = %.3lf, c = %d, size1 = %lu, size2 = %lu, sp1 = %d-%d, sp2 = %d-%d\n", 
					r, c, splices[i].size(), splices[j].size(), splices[i].front(), splices[i].back(), splices[j].front(), splices[j].back());

			if(c <= 0.50) continue;
			if(r < min_similarity) continue;

			if(local == true) vpid.push_back(PPID(PI(xi, xj), r));
			else vpid.push_back(PPID(PI(i, j), r));
		}
	}

	sort(vpid.begin(), vpid.end(), [](const PPID &x, const PPID &y){ return x.second > y.second; });

	return 0;
}

int bundle_group::test_overlap_similarity()
{
	vector<PPID> vpid;
	for(ISMI it = jindex.begin(); it != jindex.end(); it++)
	{
		const set<int> &s = it->second;
		vector<int> v(s.begin(), s.end());
		build_overlap_similarity(v, vpid, false);
	}
	return 0;
}

int bundle_group::build_overlap_similarity(const vector<int> &ss, vector<PPID> &vpid, bool local)
{
	//printf("START BUILD SIMILARITY; cfg.max_num_junctions_to_combine = %d, ss.size() =%lu\n", cfg.max_num_junctions_to_combine, ss.size());
	for(int xi = 0; xi < ss.size(); xi++)
	{
		int i = ss[xi];
		for(int xj = 0; xj < ss.size(); xj++)
		{
			int j = ss[xj];
			if(i >= j) continue;

			assert(gset[i].chrm == gset[j].chrm);
			assert(gset[i].strand == gset[j].strand);

			int32_t len = get_overlapped_length(jmaps[i], jmaps[j]);
			int32_t len1 = get_total_length(jmaps[i]);
			int32_t len2 = get_total_length(jmaps[j]);

			/*
			if(c <= 0.50) continue;
			if(r < min_similarity) continue;

			if(local == true) vpid.push_back(PPID(PI(xi, xj), r));
			else vpid.push_back(PPID(PI(i, j), r));
			*/

			double o1 = len * 100.0 / len1;
			double o2 = len * 100.0 / len2;
			double oo = o1 > o2 ? o1 : o2;

			vector<int32_t> vv(splices[i].size() + splices[j].size(), 0);
			vector<int32_t>::iterator it = set_intersection(splices[i].begin(), splices[i].end(), splices[j].begin(), splices[j].end(), vv.begin());
			int c = it - vv.begin();
			int small = splices[i].size() < splices[j].size() ? splices[i].size() : splices[j].size();
			double r = c * 1.0 / small;

			if(cfg.verbose >= 2) printf("combined-similarity: r = %.3lf, c = %d, sp1 = %lu, sp2 = %lu, len1 = %d, len2 = %d, o1 = %.1lf, o2 = %.1lf, oo = %.1lf\n", 
					r, c, splices[i].size(), splices[j].size(), len1, len2, o1, o2, oo);

			if(oo < 0.75) continue;
			if(local == true) vpid.push_back(PPID(PI(xi, xj), oo));
			else vpid.push_back(PPID(PI(i, j), oo));

		}
	}

	sort(vpid.begin(), vpid.end(), [](const PPID &x, const PPID &y){ return x.second > y.second; });

	return 0;
}

int bundle_group::augment_disjoint_set(const vector<PPID> &vpid, disjoint_set &ds)
{
	for(int i = 0; i < vpid.size(); i++)
	{
		int x = vpid[i].first.first;
		int y = vpid[i].first.second;
		int px = ds.find_set(x);
		int py = ds.find_set(y);
		if(px == py) continue;

		int sx = ds.get_size(px);
		int sy = ds.get_size(py);
		if(sx >= cfg.max_group_size) continue;
		if(sy >= cfg.max_group_size) continue;

		ds.link(px, py);

		int q = ds.find_set(px);
		assert(q == ds.find_set(py));
		ds.set_size(q, sx + sy);
	}
	return 0;
}

int bundle_group::build_groups(disjoint_set &ds)
{
	map<int, int> mm;
	for(int i = 0; i < gset.size(); i++)
	{
		int p = ds.find_set(i);
		int s = ds.get_size(p);

		if(mm.find(p) == mm.end())
		{
			vector<int> gv;
			gv.push_back(i);
			mm.insert(pair<int, int>(p, gvv.size()));
			gvv.push_back(std::move(gv));
		}
		else
		{
			int k = mm[p];
			gvv[k].push_back(i);
		}
	}
	return 0;
}

int bundle_group::filter(const set<int> &s, disjoint_set &ds, vector<int> &ss)
{
	for(auto &z: s)
	{
		if(grouped[z] == true) continue;
		int p = ds.find_set(z);
		if(ds.get_size(p) >= cfg.max_group_size) 
		{
			grouped[z] = true;
			continue;
		}
		ss.push_back(z);
	}
	return 0;
}

int bundle_group::stats(disjoint_set &ds, int r)
{
	int cnt = 0;
	for(int z = 0; z < gset.size(); z++)
	{
		int p = ds.find_set(z);
		int s = ds.get_size(p);
		if(grouped[z] == true) assert(s >= cfg.max_group_size);
		if(ds.get_size(p) >= cfg.max_group_size) cnt++;
	}

	printf("bundle group stats: round %d, chrm %s, rid %d, strand %c, total %d / %lu graphs grouped\n", r, chrm.c_str(), rid, strand, cnt, gset.size());

	map<int, int> m;
	for(int k = 0; k < gvv.size(); k++)
	{
		int n = gvv[k].size();
		if(m.find(n) == m.end()) m.insert(pair<int, int>(n, 1));
		else m[n]++;
	}

	for(map<int, int>::iterator it = m.begin(); it != m.end(); it++)
	{
		printf("bundle group stats: round %d, chrm %s, rid %d, strand %c, total %d graphs with combined %d graphs\n", r, chrm.c_str(), rid, strand, it->second, it->first);
	}
	return 0;
}

int bundle_group::print()
{
	for(int k = 0; k < gvv.size(); k++)
	{
		printf("combined-graph with %lu children: ", gvv[k].size());
		for(int i = 0; i < gvv[k].size(); i++)
		{
			int g = gvv[k][i];
			printf(" graph %d, rid = %s, splices = %lu, hits = %lu\n", g, gset[g].gid.c_str(), splices[g].size(), gset[g].hits.size());
		}
	}
	printf("\n");
	return 0;
}
