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
#include "boost/pending/disjoint_sets.hpp"

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
	build_splice_map();
	build_similarity();
	combine_graphs();
	//stats();
	mis.clear();
	vpid.clear();
	return 0;
}

int combined_group::build_splice_map()
{
	mis.clear();
	for(int k = 0; k < gset.size(); k++)
	{
		for(int i = 0; i < gset[k].splices.size(); i++)
		{
			int32_t p = gset[k].splices[i];
			MISI::iterator it = mis.find(p);
			if(it == mis.end())
			{
				set<int> s;
				s.insert(k);
				mis.insert(PISI(p, s));
			}
			else
			{
				it->second.insert(k);
			}
		}
	}
	return 0;
}

int combined_group::build_similarity()
{
	// maintain the similarity between any two graphs
	vector< set<int> > gmap(gset.size());		// success; into vpid
	vector< set<int> > fmap(gset.size());		// compared but failed

	for(MISI::iterator mi = mis.begin(); mi != mis.end(); mi++)
	{
		set<int> &ss = mi->second;
		vector<int> v(ss.begin(), ss.end());
		for(int xi = 0; xi < v.size(); xi++)
		{
			int i = v[xi];
			if(gset[i].junctions.size() > cfg.max_junctions_combine) continue;
			for(int xj = 0; xj < v.size(); xj++)
			{
				int j = v[xj];
				if(i >= j) continue;
				if(gset[j].junctions.size() > cfg.max_junctions_combine) continue;

				assert(gset[i].chrm == gset[j].chrm);
				assert(gset[i].strand == gset[j].strand);

				if(gmap[i].find(j) != gmap[i].end()) continue;
				if(fmap[i].find(j) != fmap[i].end()) continue;

				int c = gset[j].get_overlapped_splice_positions(gset[i].splices);
				double r = c * 1.0 / (gset[i].splices.size() + gset[j].splices.size() - c);

				// TODO parameters
				bool b = true;
				if(c <= 1.50) b = false;
				if(r < cfg.min_merge_similarity) b = false;

				//printf("graph-similarity: b = %c, r = %.3lf, c = %d, size1 = %lu, size2 = %lu\n", b ? 'T' : 'F', r, c, gset[i].splices.size(), gset[j].splices.size());

				if(b == true)
				{
					gmap[i].insert(j);
					vpid.push_back(PPID(PI(i, j), r));
				}
				else
				{
					fmap[i].insert(j);
				}
			}
		}
	}

	return 0;
}

int combined_group::combine_graphs()
{
	// maintain num_combined
	vector<int> csize(gset.size(), 0);
	for(int i = 0; i < gset.size(); i++)
	{
		csize[i] = gset[i].num_combined;
	}

	// disjoint map, maintain gvv
	vector<int> rank(gset.size(), -1);
	vector<int> parent(gset.size(), -1);

	disjoint_sets<int*, int*> ds(&rank[0], &parent[0]);
	for(int k = 0; k < gset.size(); k++) ds.make_set(k);

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

	map<int, int> mm;
	gvv.clear();
	for(int i = 0; i < gset.size(); i++)
	{
		int p = ds.find_set(i);
		if(mm.find(p) == mm.end())
		{
			vector<int> gv;
			gv.push_back(i);
			mm.insert(pair<int, int>(p, gvv.size()));
			gvv.push_back(gv);
		}
		else
		{
			int k = mm[p];
			gvv[k].push_back(i);
		}
	}
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

bool compare_graph_similarity(const PPID &x, const PPID &y)
{
	return x.second > y.second;
}
