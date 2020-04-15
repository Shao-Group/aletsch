#include "incubator.h"
#include "generator.h"
#include "filter.h"
#include "cluster.h"
#include "parameters.h"
#include "scallop.h"
#include "graph_reviser.h"
#include "bridge_solver.h"
#include "essential.h"
#include "constants.h"

#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <thread>
#include <ctime>
#include <boost/asio/post.hpp>
#include <boost/asio/thread_pool.hpp>
#include <boost/pending/disjoint_sets.hpp>

incubator::incubator(const parameters &c)
	: cfg(c)
{
	g2g.resize(3);
}

int incubator::resolve()
{
	time_t mytime;
	mytime = time(NULL);
	printf("\nStep 1: generate graphs for individual bam/sam files, %s\n", ctime(&mytime));
	generate();

	mytime = time(NULL);
	printf("Step 2: merge splice graphs, %s\n", ctime(&mytime));
	merge();

	mytime = time(NULL);
	printf("Step 3: assemble merged splice graphs, %s\n", ctime(&mytime));
	assemble();

	mytime = time(NULL);
	printf("Step 4: filter and output assembled transcripts, %s\n", ctime(&mytime));
	postprocess();

	return 0;
}

int incubator::generate()
{
	ifstream fin(cfg.input_bam_list.c_str());
	if(fin.fail())
	{
		printf("cannot open input-bam-list-file %s\n", cfg.input_bam_list.c_str());
		exit(0);
	}

	int num_threads = cfg.max_threads;
	if(cfg.single_sample_multiple_threading == true) num_threads = cfg.max_threads / 2;
	
	boost::asio::thread_pool pool(num_threads);			// thread pool
	mutex mylock;										// lock for trsts

	char line[102400];
	while(fin.getline(line, 10240, '\n'))
	{
		string file(line);
		if(file.size() == 0) continue;
		boost::asio::post(pool, [this, &mylock, file]{ this->generate(file, this->groups, mylock); });
	}

	pool.join();
	print_groups();

	time_t mytime = time(NULL);
	printf("finish processing all individual samples, %s\n", ctime(&mytime));

	fin.close();
	return 0;
}

int incubator::merge()
{
	boost::asio::thread_pool pool(cfg.max_threads); // thread pool

	for(int k = 0; k < groups.size(); k++)
	{
		combined_group &gp = groups[k];
		boost::asio::post(pool, [&gp]{ gp.resolve(); });
	}

	pool.join();
	print_groups();

	time_t mytime = time(NULL);
	printf("finish merging all splice graphs, %s\n", ctime(&mytime));
	return 0;	
}

int incubator::assemble()
{
	init_transcript_sets();

	boost::asio::thread_pool pool(cfg.max_threads);
	mutex mylock;

	int instance = 0;
	for(int i = 0; i < groups.size(); i++)
	{
		for(int k = 0; k < groups[i].gvv.size(); k++)
		{
			const vector<int> &v = groups[i].gvv[k];
			if(v.size() == 0) continue;
			vector<combined_graph*> gv;
			for(int j = 0; j < v.size(); j++) gv.push_back(&(groups[i].gset[v[j]]));
			boost::asio::post(pool, [this, gv, instance, &mylock]{ this->assemble(gv, instance, mylock); });
			instance++;
		}
	}
	pool.join();

	time_t mytime = time(NULL);
	printf("finish assembling all merged graphs, %s\n", ctime(&mytime));
	return 0;
}

int incubator::postprocess()
{
	ofstream fout(cfg.output_gtf_file.c_str());
	if(fout.fail())
	{
		printf("cannot open output-get-file %s\n", cfg.output_gtf_file.c_str());
		exit(0);
	}

	boost::asio::thread_pool pool(cfg.max_threads); // thread pool
	mutex mylock;									// lock for trsts

	for(int i = 0; i < tss.size(); i++)
	{
		transcript_set &ts = tss[i];
		boost::asio::post(pool, [this, &ts, &mylock, &fout]{ this->postprocess(ts, fout, mylock); });
	}

	pool.join();
	fout.close();

	time_t mytime = time(NULL);
	printf("finish filtering and reporting final assembled transcripts, %s\n", ctime(&mytime));

	return 0;
}

int incubator::generate(const string &file, vector<combined_group> &gv, mutex &mylock)
{	
	vector<combined_graph> v;
	parameters c = cfg;
	generator gt(file, v, c);
	gt.resolve();

	mylock.lock();
	for(int k = 0; k < v.size(); k++)
	{
		v[k].print(k);

		string chrm = v[k].chrm;
		char c = v[k].strand;
		int s = 0;
		if(c == '+') s = 1;
		if(c == '-') s = 2;
		if(g2g[s].find(chrm) == g2g[s].end())
		{
			combined_group gp(chrm, c, cfg);
			gp.add_graph(std::move(v[k]));
			g2g[s].insert(pair<string, int>(chrm, gv.size()));
			gv.push_back(std::move(gp));
		}
		else
		{
			gv[g2g[s][chrm]].add_graph(std::move(v[k]));
		}
	}
	mylock.unlock();

	time_t mytime = time(NULL);
	printf("finish processing individual sample %s, %s", file.c_str(), ctime(&mytime));
	return 0;
}


int incubator::assemble(vector<combined_graph*> gv, int instance, mutex &mylock)
{
	if(gv.size() == 0) return 0;

	transcript_set ts;
	int subindex = 0;

	if(gv.size() == 1)
	{
		assemble(gv, instance, subindex, ts);
		ts.increase_count(1);
		store_transcripts(ts, mylock);
		return 0;
	}

	combined_graph cx;
	resolve_cluster(gv, cx);
	cx.set_gid(instance, subindex++);
	transcript_set ts1;
	assemble(cx, ts1);

	int n = gv.size();
	set<int> ss;
	for(int i = 1; i <= 5; i++)
	{
		int k = n * i / 5;
		if(k == 1) continue;
		if(k >= n) k = n;
		ss.insert(k);
	}

	for(auto &k: ss)
	{
		vector<vector<combined_graph*>> gvv(k);
		for(int i = 0; i < gv.size(); i++)
		{
			int j = i % k;
			gvv[j].push_back(gv[i]);
		}

		transcript_set tsk;
		for(int i = 0; i < k; i++)
		{
			assemble(gvv[i], instance, subindex++, tsk);
		}

		if(k == n) tsk.add(ts1, TRANSCRIPT_COUNT_ADD_COVERAGE_NUL);
		
		ts.add(tsk, TRANSCRIPT_COUNT_MAX_COVERAGE_MAX);
	}

	store_transcripts(ts, mylock);

	for(int i = 0; i > gv.size(); i++) gv[i]->clear();

	return 0;
}

int incubator::assemble(vector<combined_graph*> gv, int instance, int subindex, transcript_set &ts)
{
	if(gv.size() <= 0) return 0;

	if(gv.size() == 1)
	{
		combined_graph &cb = *(gv[0]);
		cb.set_gid(instance, subindex);
		assemble(cb, ts);
	}
	else
	{
		combined_graph cb;
		cb.combine(gv);
		cb.copy_meta_information(*(gv[0]));
		cb.set_gid(instance, subindex);
		assemble(cb, ts);
	}
	return 0;
}

int incubator::assemble(combined_graph &cb, transcript_set &ts)
{
	// rebuild splice graph
	splice_graph gx;
	cb.build_splice_graph(gx, cfg);
	gx.build_vertex_index();

	phase_set px = cb.ps;

	map<int32_t, int32_t> smap, tmap;
	group_start_boundaries(gx, smap, cfg.max_group_boundary_distance);
	group_end_boundaries(gx, tmap, cfg.max_group_boundary_distance);
	px.project_boundaries(smap, tmap);

	refine_splice_graph(gx);

	hyper_set hx(gx, px);
	hx.filter_nodes(gx);

	/*
	cb.print(0);
	printf("==\n");
	gx.print();
	printf("==\n");
	hx.print_nodes();
	printf("==\n");
	*/

	gx.gid = cb.gid;
	scallop sx(gx, hx, cfg);
	sx.assemble();

	int z = 0;
	for(int k = 0; k < sx.trsts.size(); k++)
	{
		transcript &t = sx.trsts[k];
		if(t.exons.size() <= 1) continue;
		t.RPKM = 0;
		t.count = 1;
		z++;
		ts.add(t, TRANSCRIPT_COUNT_ADD_COVERAGE_ADD);
		//t.write(cout);
	}

	printf("assemble combined-graph %s, %d assembled transcripts: ", cb.gid.c_str(), z);
	cb.print(0);

	return 0;
}

int incubator::resolve_cluster(vector<combined_graph*> gv, combined_graph &cb)
{
	if(gv.size() <= 1) return 0;

	// construct combined graph
	cb.copy_meta_information(*(gv[0]));
	cb.combine(gv);

	// rebuild splice graph
	splice_graph gx;
	cb.build_splice_graph(gx, cfg);
	gx.build_vertex_index();

	// collect and bridge all unbridged pairs
	vector<pereads_cluster> vc;
	vector<PI> index(gv.size());
	int length_low = 999;
	int length_high = 0;
	for(int k = 0; k < gv.size(); k++)
	{
		combined_graph &gt = *(gv[k]);
		if(gt.sp.insertsize_low < length_low) length_low = gt.sp.insertsize_low;
		if(gt.sp.insertsize_high > length_high) length_high = gt.sp.insertsize_high;
		index[k].first = vc.size();
		vc.insert(vc.end(), gt.vc.begin(), gt.vc.end());
		index[k].second = vc.size();
	}

	bridge_solver br(gx, vc, cfg, length_low, length_high);
	br.build_phase_set(cb.ps);

	// resolve individual graphs
	for(int i = 0; i < gv.size(); i++)
	{
		combined_graph g1;
		for(int k = index[i].first; k < index[i].second; k++)
		{
			if(br.opt[k].type < 0) continue;
			g1.append(vc[k], br.opt[k]);
		}
		gv[i]->combine(&g1);
		gv[i]->vc.clear();
	}
	return 0;
}

int incubator::postprocess(const transcript_set &ts, ofstream &fout, mutex &mylock)
{
	vector<transcript> v = ts.get_transcripts(2);

	cluster cs(v, cfg);
	cs.solve();

	filter ft(/*v, */cs.cct, cfg);
	ft.join_single_exon_transcripts();
	ft.filter_length_coverage();

	stringstream ss;
	for(int i = 0; i < ft.trs.size(); i++)
	{
		transcript &t = ft.trs[i];
		t.write(ss);
	}

	mylock.lock();
	const string &s = ss.str();
	fout.write(s.c_str(), s.size());
	mylock.unlock();

	return 0;
}

int incubator::init_transcript_sets()
{
	tss.resize(groups.size());
	assert(g2g.size() == 3);
	string ss = ".+-";
	for(int i = 0; i < 3; i++)
	{
		for(auto &x : g2g[i])
		{
			string s = x.first;
			int k = x.second;
			assert(k >= 0 && k < tss.size());
			tss[k].chrm = s;
			tss[k].strand = ss[i];
		}
	}
	return 0;
}

int incubator::store_transcripts(const transcript_set &ts, mutex &mylock)
{
	if(ts.mt.size() == 0) return 0;

	mylock.lock();

	int k = 0;
	if(ts.strand == '+') k = 1;
	if(ts.strand == '-') k = 2;

	assert(g2g[k].find(ts.chrm) != g2g[k].end());
	int z = g2g[k][ts.chrm];

	tss[z].add(ts, 2, TRANSCRIPT_COUNT_ADD_COVERAGE_ADD);

	mylock.unlock();
	return 0;
}

int incubator::print_groups()
{
	for(int k = 0; k < groups.size(); k++)
	{
		printf("group %d (chrm = %s, strand = %c) contains %lu graphs (%lu merged graphs)\n", k, groups[k].chrm.c_str(), groups[k].strand, groups[k].gset.size(), groups[k].gvv.size());
		//groups[k].print();
	}
	return 0;
}
