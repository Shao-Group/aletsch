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

	// test the time / memory up to 2nd step
	//return 0;

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

	boost::asio::thread_pool pool(cfg.max_threads / 2); // thread pool
	mutex mylock;										// lock for trsts

	char line[102400];
	while(fin.getline(line, 10240, '\n'))
	{
		string file(line);
		if(file.size() == 0) continue;
		boost::asio::post(pool, [this, &mylock, file]{ generate_single(file, this->groups, mylock, this->g2g, this->cfg); });
		//generate_single(file, this->groups, mylock, this->g2g, this->cfg);
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
	boost::asio::thread_pool pool(cfg.max_threads); // thread pool
	mutex mylock;									// lock for trsts

	int instance = 0;
	for(int i = 0; i < groups.size(); i++)
	{
		for(int k = 0; k < groups[i].gvv.size(); k++)
		{
			const vector<int> &v = groups[i].gvv[k];
			if(v.size() == 0) continue;

			if(v.size() == 1)
			{
				combined_graph &cb = groups[i].gset[v[0]];
				boost::asio::post(pool, [this, &cb, instance, &mylock]{ assemble_single(cb, instance, this->tset, mylock, this->cfg); });
			}
			else
			{
				assert(v.size() >= 2);
				vector<combined_graph*> gv;
				for(int j = 0; j < v.size(); j++) gv.push_back(&(groups[i].gset[v[j]]));
				boost::asio::post(pool, [this, gv, instance, &mylock]{ assemble_cluster(gv, instance, this->tset, mylock, this->cfg); });
			}
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

	typedef map<size_t, vector<transcript>>::iterator MIT;
	auto &trsts = tset.mt;

	int index = 0;
	for(;;)
	{
		if(index >= trsts.size()) break;

		MIT m1 = trsts.begin();
		std::advance(m1, index);
		MIT m2 = m1;
		if(trsts.size() - index <= 1000) 
		{
			m2 = trsts.end();
			index = trsts.size();
		}
		else 
		{
			std::advance(m2, 1000);
			index += 1000;
		}

		boost::asio::post(pool, [m1, m2, &fout, &mylock, this]
				{ 
					stringstream ss;
					for(MIT x = m1; x != m2; x++)
					{
						auto &v = x->second;
						cluster cs(v, this->cfg);
						cs.solve();

						filter ft(cs.cct, this->cfg);
						ft.join_single_exon_transcripts();
						ft.filter_length_coverage();

						for(int i = 0; i < ft.trs.size(); i++)
						{
							transcript &t = ft.trs[i];
							t.RPKM = 0;
							t.write(ss);
						}
					}
					mylock.lock();
					const string &s = ss.str();
					fout.write(s.c_str(), s.size());
					ss.str("");
					mylock.unlock();
				}
		);
	}

	pool.join();
	fout.close();

	time_t mytime = time(NULL);
	printf("finish filtering and reporting final assembled transcripts, %s\n", ctime(&mytime));

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

int generate_single(const string &file, vector<combined_group> &gv, mutex &mylock, vector< map<string, int> > &m, const parameters &cfg)
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
		if(m[s].find(chrm) == m[s].end())
		{
			combined_group gp(chrm, c, cfg);
			gp.add_graph(std::move(v[k]));
			m[s].insert(pair<string, int>(chrm, gv.size()));
			gv.push_back(std::move(gp));
		}
		else
		{
			gv[m[s][chrm]].add_graph(std::move(v[k]));
		}
	}
	mylock.unlock();

	time_t mytime = time(NULL);
	printf("finish processing individual sample %s, %s", file.c_str(), ctime(&mytime));
	return 0;
}

int assemble(combined_graph &cb, transcript_set &vt, const parameters &cfg, bool group_boundary)
{
	// rebuild splice graph
	splice_graph gx;
	cb.build_splice_graph(gx);
	gx.build_vertex_index();

	phase_set px = cb.ps;
	if(group_boundary == true)
	{
		map<int32_t, int32_t> smap, tmap;
		group_start_boundaries(gx, smap, cfg.max_group_boundary_distance);
		group_end_boundaries(gx, tmap, cfg.max_group_boundary_distance);
		px.project_boundaries(smap, tmap);
	}

	refine_splice_graph(gx);

	hyper_set hx(gx, px);
	hx.filter_nodes(gx);

	// assemble 
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
		vt.add(t, ADD_TRANSCRIPT_COVERAGE_SUM);
	}

	printf("assemble combined-graph %s, %d assembled transcripts: ", cb.gid.c_str(), z);
	cb.print(0);

	return 0;
}

int assemble(vector<combined_graph*> gv, int instance, int subindex, transcript_set &vt, const parameters &cfg)
{
	if(gv.size() <= 0) return 0;

	if(gv.size() == 1)
	{
		combined_graph &cb = *(gv[0]);
		cb.set_gid(instance, subindex);
		assemble(cb, vt, cfg, false);
	}
	else
	{
		combined_graph cb;
		cb.combine(gv);
		cb.copy_meta_information(*(gv[0]));
		cb.set_gid(instance, subindex);
		assemble(cb, vt, cfg, true);
	}
	return 0;
}

int assemble_single(combined_graph &cb, int instance, transcript_set &ts, mutex &mylock, const parameters &cfg)
{
	cb.set_gid(instance, 0);
	transcript_set vt;
	assemble(cb, vt, cfg, false);

	mylock.lock();
	ts.add(vt, ADD_TRANSCRIPT_COVERAGE_SUM);
	mylock.unlock();

	cb.clear();
	return 0;
}

int assemble_cluster(vector<combined_graph*> gv, int instance, transcript_set &ts, mutex &mylock, const parameters &cfg)
{
	assert(gv.size() >= 2);
	transcript_set tts;

	int subindex = 0;
	combined_graph cx;
	resolve_cluster(gv, cx, cfg);
	cx.set_gid(instance, subindex++);
	transcript_set vt0;
	assemble(cx, vt0, cfg, true);

	// sample d points from [1, n]
	set<int> ss;
	/*
	int d = 4;
	int n = (1 + gv.size()) / 2 - 1;
	if(d > n) d = n;
	ss.insert(1);
	ss.insert(n + 1);
	for(int i = 1; i < d; i++) ss.insert(1 + i * n / d);
	*/

	ss.insert(1);

	for(auto &k: ss)
	{
		//if(k < 1 || k > n + 1) continue;
		transcript_set vt;
		for(int i = 0; i <= gv.size() / k; i++)
		{
			vector<combined_graph*> gv1;
			for(int j = 0; j < k; j++)
			{
				if(i * k + j < gv.size()) gv1.push_back(gv[i * k + j]);
			}
			assemble(gv1, instance, subindex++, vt, cfg);
		}
		if(k == 1) vt.add(vt0, ADD_TRANSCRIPT_COVERAGE_NAN);
		tts.add(vt, 2, ADD_TRANSCRIPT_COVERAGE_NAN);
	}

	//printf("instance = %d, %lu predicted transcripts\n", instance, tts.get_transcripts(1).size());

	mylock.lock();
	//tts.print();
	ts.add(tts, ADD_TRANSCRIPT_COVERAGE_SUM);
	mylock.unlock();

	for(int i = 0; i > gv.size(); i++) gv[i]->clear();
	return 0;
}

int resolve_cluster(vector<combined_graph*> gv, const parameters &cfg)
{
	combined_graph cb;
	resolve_cluster(gv, cb, cfg);
	return 0;
}

int resolve_cluster(vector<combined_graph*> gv, combined_graph &cb, const parameters &cfg)
{
	if(gv.size() <= 1) return 0;

	// construct combined graph
	cb.copy_meta_information(*(gv[0]));
	cb.combine(gv);

	// rebuild splice graph
	splice_graph gx;
	cb.build_splice_graph(gx);
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
