#include "combine.h"
#include "generate.h"
#include "undirected_graph.h"
#include "boost/pending/disjoint_sets.hpp"
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <thread>
#include <boost/asio/thread_pool.hpp>
#include <boost/asio/post.hpp>

incubator::incubator(int m, int t)
{
	max_threads = t;
	max_combined = m;
	g2g.resize(3);
}

int incubator::load(const string &file)
{
	ifstream fin(file.c_str());
	if(fin.fail())
	{
		printf("cannot open file %s\n", file.c_str());
		exit(0);
	}

	vector< vector<string> > files(max_threads);

	char line[102400];
	int index = 0;
	while(fin.getline(line, 10240, '\n'))
	{
		string s(line);
		if(s.size() == 0) continue;
		int k = index % max_threads;
		files[k].push_back(s);
		index++;
	}

	mutex mylock;								// lock for trsts
	vector<thread> threads;
	for(int k = 0; k < files.size(); k++)
	{
		printf("thread %d processes %lu files\n", k, files[k].size());
		threads.emplace_back(load_multiple, files[k], std::ref(groups), std::ref(mylock), std::ref(g2g));
	}
	for(int k = 0; k < threads.size(); k++)
	{
		threads[k].join();
	}

	print_groups();
	return 0;
}

int incubator::merge(double ratio)
{
	boost::asio::thread_pool pool(max_threads); // thread pool

	for(int k = 0; k < groups.size(); k++)
	{
		combined_group &gp = groups[k];
		boost::asio::post(pool, [&gp, this, ratio]{ gp.resolve(this->max_combined, ratio); });
	}

	pool.join();
	print_groups();

	return 0;	
}

int incubator::binary_merge(const string &file)
{
	ifstream fin(file.c_str());
	if(fin.fail())
	{
		printf("cannot open file %s\n", file.c_str());
		exit(0);
	}

	vector<string> files;

	char line[102400];
	while(fin.getline(line, 10240, '\n'))
	{
		string s(line);
		if(s.size() == 0) continue;
		files.push_back(s);
	}

	vector<combined_graph> vc;
	binary_merge(files, 0, files.size(), vc, true);

	fixed.insert(fixed.end(), vc.begin(), vc.end());
	return 0;
}

int incubator::binary_merge(const vector<string> &files, int low, int high, vector<combined_graph> &vc, bool last)
{
	//printf("binary merge from %d to %d (total = %lu)\n", low, high, files.size());
	vc.clear();

	if(low >= high) return 0;

	if(low + 1 == high)
	{
		string file = files[low];
		load_single(file, vc);
		printf("create %lu combined-graphs for file %d (%s)\n", vc.size(), low, files[low].c_str());
		return 0;
	}

	int mid = (low + high) / 2;

	vector<combined_graph> vc1;
	vector<combined_graph> vc2;
	binary_merge(files, low, mid, vc1, false);
	binary_merge(files, mid, high, vc2, false);

	int n1 = vc1.size();
	int n2 = vc2.size();

	vc1.insert(vc1.end(), vc2.begin(), vc2.end());

	merge(vc1, vc, last);

	printf("merge final with %lu <- %d + %d (fixed = %lu) combined-graphs for files [%d, %d]\n", vc.size(), n1, n2, fixed.size(), low, high - 1);

	if(mdir != "" && high - low >= 10)
	{
		char file[10240];
		sprintf(file, "%s/graph-%d-%d.gr", mdir.c_str(), low, high - 1);
		write(file, true);
	}

	return 0;
}

int incubator::merge(const vector<combined_graph> &grset, vector<combined_graph> &vc, bool last)
{
	// maintain num_combined
	vector<int> csize(grset.size(), 0);
	for(int i = 0; i < grset.size(); i++)
	{
		csize[i] = grset[i].num_combined;
	}

	// disjoint map, maintain clusters
	vector<int> rank(grset.size(), -1);
	vector<int> parent(grset.size(), -1);

	disjoint_sets<int*, int*> ds(&rank[0], &parent[0]);
	for(int k = 0; k < grset.size(); k++)
	{
		ds.make_set(k);
	}

	// maintain the similarity between any two graphs
	vector< map<int, double> > gmap(grset.size());

	vector<PID> vpid;

	MISI mis;
	build_splice_map(grset, mis);
	for(MISI::iterator mi = mis.begin(); mi != mis.end(); mi++)
	{
		set<int> &ss = mi->second;
		vector<int> v(ss.begin(), ss.end());
		for(int xi = 0; xi < v.size(); xi++)
		{
			int i = v[xi];
			assert(grset[i].num_combined < max_combined);
			for(int xj = 0; xj < v.size(); xj++)
			{
				int j = v[xj];
				if(i >= j) continue;

				assert(grset[j].num_combined < max_combined);

				if(grset[i].chrm != grset[j].chrm) continue;
				if(grset[i].strand != grset[j].strand) continue;

				if(gmap[i].find(j) != gmap[i].end()) continue;

				int c = grset[j].get_overlapped_splice_positions(grset[i].splices);
				if(c <= 1.5) continue;

				/*
				double r1 = c * 1.0 / grset[i].splices.size();
				double r2 = c * 1.0 / grset[j].splices.size();
				double rx = r1 < r2 ? r1 : r2;
				double ry = r1 > r2 ? r1 : r2;
				double r = 2 * rx + 0.5 * ry;
				*/

				// TODO parameter
				double r = 2.0 * c / (grset[i].splices.size() + grset[j].splices.size());
				if(last == false && r < 0.8) continue;

				//printf("r1 = %.3lf, r2 = %.3lf, r = %.3lf, size1 = %lu, size2 = %lu\n", r1, r2, r, grset[i].splices.size(), grset[j].splices.size());

				gmap[i].insert(pair<int, double>(j, r));
				vpid.push_back(PID(PI(i, j), r));
			}
		}
	}

	sort(vpid.begin(), vpid.end(), compare_graph_overlap);

	for(int i = 0; i < vpid.size(); i++)
	{
		int x = vpid[i].first.first;
		int y = vpid[i].first.second;
		double r = vpid[i].second;
		assert(x < y);

		int px = ds.find_set(x);
		int py = ds.find_set(y);

		if(px == py) continue;
		if(csize[px] >= max_combined) continue;
		if(csize[py] >= max_combined) continue;

		int sum = csize[px] + csize[py]; 

		printf("combine graph %d (#splices = %lu) and %d (#splices = %lu) with score = %.3lf: %d + %d -> %d\n", 
				x, grset[x].splices.size(), y, grset[y].splices.size(), r, csize[px], csize[py], sum);

		ds.link(px, py);
		int q = ds.find_set(px);
		assert(q == ds.find_set(py));
		assert(q == px || q == py);
		csize[q] = sum;
	}

	vector<combined_graph> cc(grset.size());
	vector<int> pp(grset.size(), -1);
	for(int i = 0; i < grset.size(); i++)
	{
		int p = ds.find_set(i);
		pp[i] = p;
		if(p == i) cc[i] = grset[i];
	}

	for(int i = 0; i < grset.size(); i++)
	{
		if(pp[i] == i) continue;
		int p = pp[i];
		cc[p].combine(grset[i]);
	}

	for(int i = 0; i < grset.size(); i++)
	{
		if(cc[i].num_combined >= max_combined)
		{
			assert(pp[i] == i);
			fixed.push_back(cc[i]);
		}
		else if(i == pp[i])
		{
			vc.push_back(cc[i]);
		}
	}
	return 0;
}

int incubator::write(const string &file, bool headers)
{
	ofstream fout(file.c_str());
	if(fout.fail()) exit(1);

	int os = 0;
	vector<int> offset;
	for(int k = 0; k < groups.size(); k++)
	{
		offset.push_back(os);
		os += groups[k].mset.size();
	}

	mutex mylock;								// lock for trsts
	boost::asio::thread_pool pool(max_threads); // thread pool
	for(int k = 0; k < groups.size(); k++)
	{
		combined_group &gp = groups[k];
		int os = offset[k];
		boost::asio::post(pool, [&gp, &fout, &mylock, os]{ gp.write(mylock, fout, os); });
	}

	pool.join();


	/*
	int index = 0;
	for(int k = 0; k < groups.size(); k++)
	{
		for(int j = 0; j < groups[k].mset.size(); j++)
		{
			groups[k].mset[j].write(fout, index++, headers);
		}
	}
	*/

	fout.flush();
	fout.close();
	return 0;
}

int incubator::build_splice_map(const vector<combined_graph> &grset, MISI &mis)
{
	mis.clear();
	for(int k = 0; k < grset.size(); k++)
	{
		for(int i = 0; i < grset[k].splices.size(); i++)
		{
			int32_t p = grset[k].splices[i];
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

int incubator::print()
{
	for(int k = 0; k < fixed.size(); k++) fixed[k].print(k);
	return 0;
}

int incubator::print_groups()
{
	for(int k = 0; k < groups.size(); k++)
	{
		printf("group %d (chrm = %s, strand = %c) contains %lu graphs (%lu merged graphs)\n", k, groups[k].chrm.c_str(), groups[k].strand, groups[k].gset.size(), groups[k].mset.size());
	}
	return 0;
}

int incubator::analyze(const string &file)
{
	vector<combined_graph> vc;
	load_single(file, vc);
	for(int k = 0; k < vc.size(); k++)
	{
		vc[k].analyze(k);
	}
	return 0;
}

int load_multiple(const vector<string> &files, vector<combined_group> &gv, mutex &mylock, vector< map<string, int> > &m)
{	
	vector<combined_graph> v;
	for(int k = 0; k < files.size(); k++) 
	{
		//printf("load file %s\n", files[k].c_str());
		//load_single(files[k], v);
		generator gt(files[k], "", v);
		gt.resolve();
	}

	mylock.lock();
	for(int k = 0; k < v.size(); k++)
	{
		string chrm = v[k].chrm;
		char c = v[k].strand;
		int s = 0;
		if(c == '+') s = 1;
		if(c == '-') s = 2;
		if(m[s].find(chrm) == m[s].end())
		{
			combined_group gp(chrm, c);
			gp.add_graph(v[k]);
			m[s].insert(pair<string, int>(chrm, gv.size()));
			gv.push_back(gp);
		}
		else
		{
			gv[m[s][chrm]].add_graph(v[k]);
		}
	}
	mylock.unlock();
	return 0;
}

int load_single(const string &file, vector<combined_graph> &vc)
{
	ifstream fin(file.c_str());
	if(fin.fail())
	{
		printf("could not load file %s\n", file.c_str());
		exit(0);
	}

	char line[10240];
	char gid[10240];
	char chrm[10240];
	char tmp[1024];
	char strand[1024];
	int nodes;

	while(fin.getline(line, 10240, '\n'))
	{
		if(line[0] != '#') continue;
		stringstream sstr(line);
		sstr >> tmp >> gid >> chrm >> strand;

		combined_graph gr(line);
		gr.build(fin, chrm, strand[0]);
		vc.push_back(gr);
	}

	printf("loaded graphs in file %s\n", file.c_str());

	fin.close();
	return 0;
}

bool compare_graph_overlap(const PID &x, const PID &y)
{
	return x.second > y.second;
}
