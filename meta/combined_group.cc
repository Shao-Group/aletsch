/*
Part of meta-scallop
(c) 2020 by Mingfu Shao, The Pennsylvania State University
See LICENSE for licensing.
*/

#include <algorithm>
#include <sstream>
#include <fstream>
#include "combined_group.h"
#include "parameters.h"

mutex combined_group::gmutex;

combined_group::combined_group(string c, char s, const parameters &f)
	: cfg(f)
{
	chrm = c;
	strand = s;
}

int combined_group::add_graph(const combined_graph &gr)
{
	assert(gr.chrm == chrm);
	gset.push_back(gr);
	return 0;
}

int combined_group::resolve()
{
	build_splice_index();
	build_disjoint_sets();

	combine_graphs();
	//stats();
	sindex.clear();
	vpid.clear();
	return 0;
}

int combined_group::build_splice_index()
{
	sindex.clear();
	bmap.assign(gset.size(), false);
	for(int k = 0; k < gset.size(); k++)
	{
		for(int i = 0; i < gset[k].splices.size(); i++)
		{
			int32_t p = gset[k].splices[i];
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

	//for(auto &x : sindex) printf("splice map %d => %lu elements\n", x.first, x.second.size());

	return 0;
}

int combined_group::init_disjoint_sets()
{
	int n = divisors + 1;
	ranks.resize(n, NULL);
	parents.resize(n, NULL);
	dsss.resize(n, NULL);
	for(int i = 0; i < n; i++)
	{
		ranks[i] = new int[gsets.size()];
		parents[i] = new int[gsets.size()];
		dsss[i] = new dss(ranks[i], parents[i]);

		for(int k = 0; k < gsets.size(); k++)
		{
			ranks[i][k] = -1;
			parents[i][k] = -1;
			dsss[i].make_set(k);
		}
	}
	return 0;
}

int combined_group::build_disjoint_sets()
{
	boost::asio::thread_pool pool(cfg.max_threads);
	for(auto &z: sindex)
	{
		const set<int> &s = z.second;
		vector<int> v(s.begin(), s.end());
		boost::asio::post(pool, [&v]{ process_subset(v); });
	}
	pool.join();
	return 0;
}

int combined_group::process_subset(const vector<int> &ss)
{
	vector<PPID> vpid;
	build_similarity(ss, vpid);
	build_clusters(ss, vpid);
	return 0;
}

int combined_group::build_similarity(const vector<int> &ss, vector<PPID> &vpid)
{
	for(int xi = 0; xi < ss.size(); xi++)
	{
		int i = ss[xi];
		if(gset[i].junctions.size() > cfg.max_junctions_combine) continue;
		for(int xj = 0; xj < ss.size(); xj++)
		{
			int j = ss[xj];
			if(i >= j) continue;
			if(gset[j].junctions.size() > cfg.max_junctions_combine) continue;

			assert(gset[i].chrm == gset[j].chrm);
			assert(gset[i].strand == gset[j].strand);

			int c = gset[j].get_overlapped_splice_positions(gset[i].splices);
			double r = c * 1.0 / (gset[i].splices.size() + gset[j].splices.size() - c);

			//printf("graph-similarity: b = %c, r = %.3lf, c = %d, size1 = %lu, size2 = %lu\n", b ? 'T' : 'F', r, c, gset[i].splices.size(), gset[j].splices.size());

			if(c <= 1.50) continue;
			if(r < cfg.min_merge_similarity) continue;

			vpid.push_back(PPID(PI(i, j), r));
		}
	}

	sort(vpid.begin(), vpid.end(), [](const PPID &x, const PPID &y){ return x.second > y.second; });

	return 0;
}

int combined_group::build_clusters(const vector<PPID> &vpid)
{
	for(int i = 0; i < vpid.size(); i++)
	{
		int x = vpid[i].first.first;
		int y = vpid[i].first.second;
		double r = vpid[i].second;

		int n = floor(r * divisors);
		for(int k = 0; k <= n; k++)
		{


			int px = ds.find_set(x);
			int py = ds.find_set(y);

			if(px == py) continue;
			if(csize[px] >= cfg.max_merge_cluster) continue;
			if(csize[py] >= cfg.max_merge_cluster) continue;

			int sum = csize[px] + csize[py]; 

			//printf("combine graph %d (#splices = %lu) and %d (#splices = %lu) with score = %.3lf: %d + %d -> %d\n", x, gset[x].splices.size(), y, gset[y].splices.size(), r, csize[px], csize[py], sum);

			ds.link(px, py);
			int q = ds.find_set(px);
			assert(q == ds.find_set(py));
			assert(q == px || q == py);
			csize[q] = sum;
		}
	}

	// LOCK, access and write to gvv
	map<int, int> mm;
	for(int i = 0; i < ss.size(); i++)
	{
		int p = ds.find_set(i);
		if(csize[p] < cfg.max_merge_cluster) continue;

		bmap[ss[i]] = true;
		if(mm.find(p) == mm.end())
		{
			vector<int> gv;
			gv.push_back(ss[i]);
			mm.insert(pair<int, int>(p, gvv.size()));
			gvv.push_back(gv);
		}
		else
		{
			int k = mm[p];
			gvv[k].push_back(ss[i]);
		}
	}
	// UNLOCK
	return 0;
}

int combined_group::build_clusters(const vector<int> &ss, const vector<PPID> &vpid)
{
	// maintain num_combined
	vector<int> csize(ss.size(), 0);
	for(int i = 0; i < ss.size(); i++)
	{
		csize[i] = gset[ss[i]].num_combined;
	}

	// disjoint map, maintain gvv
	vector<int> rank(ss.size(), -1);
	vector<int> parent(ss.size(), -1);

	disjoint_sets<int*, int*> ds(&rank[0], &parent[0]);
	for(int k = 0; k < ss.size(); k++) ds.make_set(k);

	sort(vpid.begin(), vpid.end(), compare_graph_similarity);

	for(int i = 0; i < vpid.size(); i++)
	{
		int x = vpid[i].first.first;
		int y = vpid[i].first.second;
		double r = vpid[i].second;
		assert(x < y);

		int px = ds.find_set(x);
		int py = ds.find_set(y);

		if(px == py) continue;
		if(csize[px] >= cfg.max_merge_cluster) continue;
		if(csize[py] >= cfg.max_merge_cluster) continue;

		int sum = csize[px] + csize[py]; 

		//printf("combine graph %d (#splices = %lu) and %d (#splices = %lu) with score = %.3lf: %d + %d -> %d\n", x, gset[x].splices.size(), y, gset[y].splices.size(), r, csize[px], csize[py], sum);

		ds.link(px, py);
		int q = ds.find_set(px);
		assert(q == ds.find_set(py));
		assert(q == px || q == py);
		csize[q] = sum;
	}

	// LOCK, access and write to gvv
	map<int, int> mm;
	for(int i = 0; i < ss.size(); i++)
	{
		int p = ds.find_set(i);
		if(csize[p] < cfg.max_merge_cluster) continue;

		bmap[ss[i]] = true;
		if(mm.find(p) == mm.end())
		{
			vector<int> gv;
			gv.push_back(ss[i]);
			mm.insert(pair<int, int>(p, gvv.size()));
			gvv.push_back(gv);
		}
		else
		{
			int k = mm[p];
			gvv[k].push_back(ss[i]);
		}
	}
	// UNLOCK
	return 0;
}

int combined_group::stats()
{
	map<int, int> m;
	for(int k = 0; k < gvv.size(); k++)
	{
		int n = gvv[k].size();
		if(m.find(n) == m.end()) m.insert(pair<int, int>(n, 1));
		else m[n]++;
	}

	for(map<int, int>::iterator it = m.begin(); it != m.end(); it++)
	{
		printf("chrm %s, strand %c, total %d graphs with combined %d graphs\n", chrm.c_str(), strand, it->second, it->first);
	}
	return 0;
}

int combined_group::print()
{
	for(int k = 0; k < gvv.size(); k++)
	{
		printf("combined-graph with %lu children: ", gvv[k].size());
		for(int i = 0; i < gvv[k].size(); i++)
		{
			int g = gvv[k][i];
			printf("%s, ", gset[g].gid.c_str());
		}
		printf("\n");
	}
	printf("\n");
	return 0;
}
