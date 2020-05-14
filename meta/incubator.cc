/*
Part of meta-scallop
(c) 2020 by Mingfu Shao, The Pennsylvania State University
See LICENSE for licensing.
*/

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

incubator::~incubator()
{
	if(cfg.output_bridged_bam_dir == "") return;
	for(int i = 0; i < samples.size(); i++)
	{
		samples[i].close_bridged_bam();
	}
}

int incubator::resolve()
{
	time_t mytime;
	read_bam_list();

	int a = 0;
	while(true)
	{
		int b = a + cfg.meta_batch_size;

		if(b >= samples.size()) b = samples.size();
		if(a >= b) break;

		mytime = time(NULL);
		printf("START PROCESSING BATCH WITH %d SAMPLES, %s\n", b - a, ctime(&mytime));

		mytime = time(NULL);
		printf("step 1: generate graphs for individual bam/sam files, %s\n", ctime(&mytime));
		generate(a, b);

		mytime = time(NULL);
		printf("step 2: merge splice graphs, %s\n", ctime(&mytime));
		merge();

		mytime = time(NULL);
		printf("step 3: assemble merged splice graphs, %s\n", ctime(&mytime));
		assemble();

		a = b;
		clear();

		mytime = time(NULL);
		printf("FINISH PROCESSING BATCH\n");
	}

	mytime = time(NULL);
	printf("\nFINAL: filter and output assembled transcripts, %s\n", ctime(&mytime));

	postprocess();
	write();

	return 0;
}

int incubator::read_bam_list()
{
	ifstream fin(cfg.input_bam_list.c_str());
	if(fin.fail())
	{
		printf("cannot open input-bam-list-file %s\n", cfg.input_bam_list.c_str());
		exit(0);
	}

	char line[102400];
	while(fin.getline(line, 102400, '\n'))
	{
		string file(line);
		if(file.size() == 0) continue;
		sample_profile sp;
		sp.file_name = file;
		if(cfg.output_bridged_bam_dir != "")
		{
			sp.open_bridged_bam(cfg.output_bridged_bam_dir);
		}
		samples.push_back(sp);
	}

	for(int i = 0; i < samples.size(); i++)
	{
		samples[i].sample_id = i;
	}

	fin.close();
	return 0;
}

int incubator::clear()
{
	for(int i = 0; i < 3; i++) g2g[i].clear();
	groups.clear();
	return 0;
}

int incubator::generate(int a, int b)
{
	if(a >= b) return 0;

	int num_threads = cfg.max_threads;
	if(cfg.single_sample_multiple_threading == true) num_threads = cfg.max_threads / 2;
	if(num_threads > b - a) num_threads = b - a;
	
	boost::asio::thread_pool pool(num_threads);			// thread pool
	mutex mylock;										// lock for tsets

	for(int i = a; i < b; i++)
	{
		sample_profile &sp = samples[i];
		boost::asio::post(pool, [this, &mylock, &sp]{ this->generate(sp, this->groups, mylock); });
	}

	pool.join();
	print_groups();

	time_t mytime = time(NULL);
	printf("finish processing all individual samples, %s\n", ctime(&mytime));

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
	mutex mylock;									// lock for tsets

	strsts.resize(samples.size());
	for(int i = 0; i < tsets.size(); i++)
	{
		transcript_set &ts = tsets[i];
		boost::asio::post(pool, [this, &ts, &mylock, &fout]{ this->postprocess(ts, fout, mylock); });
	}

	pool.join();
	fout.close();

	time_t mytime = time(NULL);
	printf("finish filtering and reporting final assembled transcripts, %s\n", ctime(&mytime));

	return 0;
}

int incubator::write()
{
	if(cfg.output_gtf_dir == "") return 0;

	boost::asio::thread_pool pool(cfg.max_threads); // thread pool
	mutex mylock;									// lock for tsets

	for(int i = 0; i < strsts.size(); i++)
	{
		boost::asio::post(pool, [this, i]{ this->write(i); });
	}

	pool.join();

	time_t mytime = time(NULL);
	printf("finish writing individual assembled transcripts, %s\n", ctime(&mytime));
	return 0;
}

int incubator::generate(sample_profile &sp, vector<combined_group> &gv, mutex &mylock)
{	
	vector<combined_graph> v;
	parameters c = cfg;
	generator gt(sp, v, c);
	gt.resolve();

	mylock.lock();
	for(int k = 0; k < v.size(); k++)
	{
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
	printf("finish processing individual sample %s, %s", sp.file_name.c_str(), ctime(&mytime));
	return 0;
}

int incubator::assemble(vector<combined_graph*> gv, int instance, mutex &mylock)
{
	if(gv.size() == 0) return 0;

	transcript_set ts;
	int subindex = 0;

	if(gv.size() == 1)
	{
		gv[0]->set_gid(instance, subindex++);
		assemble(*gv[0], ts, TRANSCRIPT_COUNT_ADD_COVERAGE_ADD);
		ts.increase_count(1);
	}
	else
	{
		combined_graph cx;
		resolve_cluster(gv, cx);

		for(int i = 0; i < gv.size(); i++)
		{
			gv[i]->set_gid(instance, subindex++);
			assemble(*gv[i], ts, TRANSCRIPT_COUNT_ADD_COVERAGE_ADD);
		}

		cx.set_gid(instance, subindex++);
		assemble(cx, ts, TRANSCRIPT_COUNT_ADD_COVERAGE_NUL);
	}

	store_transcripts(ts, mylock);
	for(int i = 0; i > gv.size(); i++) gv[i]->clear();

	return 0;
}

int incubator::assemble(combined_graph &cb, transcript_set &ts, int mode)
{
	// rebuild splice graph
	splice_graph gx;
	cb.build_splice_graph(gx, cfg);
	gx.build_vertex_index();
	gx.extend_strands();

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

	/*
	gx.print();
	if(gx.num_vertices() <= 40) 
	{
		string texfile = "tex/" + gx.gid + ".tex";
		gx.draw(texfile);
	}
	*/

	scallop sx(gx, hx, cfg);
	sx.assemble();

	int z = 0;
	for(int k = 0; k < sx.trsts.size(); k++)
	{
		transcript &t = sx.trsts[k];
		if(t.exons.size() <= 1) continue;
		t.RPKM = 0;
		z++;
		ts.add(t, 1, cb.sid, mode);
		//t.write(cout);
	}

	//printf("assemble %s: %d transcripts, ", cb.gid.c_str(), z);
	//cb.print(0);

	return 0;
}

int incubator::resolve_cluster(vector<combined_graph*> gv, combined_graph &cb)
{
	assert(gv.size() >= 2);

	// construct combined graph
	cb.copy_meta_information(*(gv[0]));
	cb.combine(gv);
	cb.sid = -1;

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
		assert(gt.sid >= 0 && gt.sid < samples.size());
		sample_profile &sp = samples[gt.sid];
		if(sp.insertsize_low < length_low) length_low = sp.insertsize_low;
		if(sp.insertsize_high > length_high) length_high = sp.insertsize_high;
		index[k].first = vc.size();
		vc.insert(vc.end(), gt.vc.begin(), gt.vc.end());
		index[k].second = vc.size();
	}

	bridge_solver br(gx, vc, cfg, length_low, length_high);
	br.build_phase_set(cb.ps);

	printf("cluster-bridge, combined = %lu, ", gv.size());
	br.print();

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
	}

	// write bridged and unbridged reads
	for(int i = 0; i < gv.size(); i++)
	{
		combined_graph &gt = *(gv[i]);
		sample_profile &sp = samples[gt.sid];
		for(int k = index[i].first; k < index[i].second; k++)
		{
			if(br.opt[k].type < 0) 
			{
				write_unbridged_pereads_cluster(sp.bridged_bam, vc[k]);
			}
			else
			{
				write_bridged_pereads_cluster(sp.bridged_bam, vc[k], br.opt[k].whole);
			}
		}
	}

	// clear to release memory
	for(int i = 0; i < gv.size(); i++) gv[i]->vc.clear();

	return 0;
}

int incubator::postprocess(const transcript_set &ts, ofstream &fout, mutex &mylock)
{
	vector<transcript> v = ts.get_transcripts(2);

	// warning: ts contains mixed strands
	//cluster cs(v, cfg);
	//cs.solve();

	filter ft(v, /*cs.cct,*/ cfg);
	ft.join_single_exon_transcripts();
	ft.filter_length_coverage();

	stringstream ss;
	vector<vector<transcript>> vv(samples.size());
	for(int i = 0; i < ft.trs.size(); i++)
	{
		transcript &t = ft.trs[i];
		t.write(ss);
		pair<bool, trans_item> p = ts.query(t);
		if(p.first == false) continue;

		for(auto &k : p.second.samples)
		{
			if(k < 0) continue;
			assert(k >= 0 && k < vv.size());
			vv[k].push_back(t);
		}
	}

	mylock.lock();

	const string &s = ss.str();
	fout.write(s.c_str(), s.size());

	for(int i = 0; i < vv.size(); i++)
	{
		strsts[i].insert(strsts[i].end(), vv[i].begin(), vv[i].end());
	}

	mylock.unlock();

	return 0;
}

int incubator::write(int id)
{
	if(id < 0) return 0;
	assert(id >= 0 && id < samples.size());
	assert(samples.size() == strsts.size());

	vector<transcript> &v = strsts[id];

	stringstream ss;
	for(int i = 0; i < v.size(); i++)
	{
		transcript &t = v[i];
		t.write(ss);
	}

	size_t p = samples[id].file_name.find_last_of("/\\");
	string bfile = samples[id].file_name.substr(p + 1);

	char file[10240];
	sprintf(file, "%s/%d.gtf", cfg.output_gtf_dir.c_str(), id);

	ofstream fout(file);
	if(fout.fail()) return 0;

	const string &s = ss.str();
	fout.write(s.c_str(), s.size());

	fout.close();
	return 0;
}

int incubator::store_transcripts(const transcript_set &ts, mutex &mylock)
{
	if(ts.mt.size() == 0) return 0;

	mylock.lock();

	bool found = false;
	for(int i = 0; i < tsets.size(); i++)
	{
		if(tsets[i].chrm != ts.chrm) continue;
		tsets[i].add(ts, 2, TRANSCRIPT_COUNT_ADD_COVERAGE_ADD);
		found = true;
		break;
	}

	if(found == false)
	{
		tsets.resize(tsets.size() + 1);
		tsets[tsets.size() - 1].add(ts, 2, TRANSCRIPT_COUNT_ADD_COVERAGE_ADD);
	}

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
