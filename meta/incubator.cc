#include "incubator.h"
#include "generator.h"
#include "filter.h"
#include "cluster.h"
#include "parameters.h"
#include "scallop.h"
#include "graph_reviser.h"
#include "bridge_solver.h"
#include "essential.h"

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
				boost::asio::post(pool, [this, &cb, instance, &mylock]{ assemble_single(cb, instance, this->trsts, mylock, this->cfg); });
			}
			else
			{
				assert(v.size() >= 2);
				vector<combined_graph*> gv;
				for(int j = 0; j < v.size(); j++) gv.push_back(&(groups[i].gset[v[j]]));
				assert(gv.size() >= 2);
				boost::asio::post(pool, [this, gv, instance, &mylock]{ assemble_cluster(gv, instance, this->trsts, mylock, this->cfg); });
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

	typedef map<size_t, vector<transcript> >::iterator MIT;
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
						vector<transcript> &v = x->second;
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

int assemble_single(combined_graph &cb, int instance, map< size_t, vector<transcript> > &trsts, mutex &mylock, const parameters &cfg1)
{
	parameters cfg = cfg1;
	vector<transcript> vt;

	// setting up names
	char name[10240];
	sprintf(name, "instance.%d.0", instance);
	cb.gid = name;

	// rebuild splice graph
	splice_graph gx;
	cb.build_splice_graph(gx);

	// refine splice graph
	refine_splice_graph(gx);
	//keep_surviving_edges(gx, cfg.min_splicing_count);

	// construct hyper-set
	hyper_set hx(gx, cb.ps);
	hx.filter_nodes(gx);

	// assemble 
	gx.gid = cb.gid;
	cfg.algo = "single";
	scallop sx(gx, hx, &cfg);
	sx.assemble();

	for(int k = 0; k < sx.trsts.size(); k++)
	{
		transcript &t = sx.trsts[k];
		t.RPKM = 0;
		if(t.exons.size() <= 1) continue;
		vt.push_back(t);
	}

	printf("assemble combined-graph %s, 0 children, %lu assembled transcripts\n", cb.gid.c_str(), vt.size());

	if(vt.size() <= 0) return 0;

	mylock.lock();
	for(int k = 0; k < vt.size(); k++)
	{
		index_transcript(trsts, vt[k]);
	}
	mylock.unlock();
	return 0;
}

int assemble_cluster(vector<combined_graph*> gv, int instance, int subindex, vector<transcript> &vt, const parameters &cfg1)
{
	assert(gv.size() >= 2);

	parameters cfg = cfg1;

	// assembled transcripts and index
	map< size_t, vector<transcript> > mt;

	// construct combined graph
	combined_graph cb;
	cb.combine(gv);

	// setting up names
	char name[10240];
	sprintf(name, "instance.%d.%d.0", instance, subindex);
	cb.gid = name;
	for(int j = 0; j < gv.size(); j++)
	{
		sprintf(name, "instance.%d.%d.%d", instance, subindex, j + 1);
		gv[j]->gid = name;
	}

	//set<int32_t> rs = cb.get_reliable_splices(cfg.min_supporting_samples, 99999);

	// rebuild splice graph
	splice_graph gx;
	cb.build_splice_graph(gx);
	phase_set px = cb.ps;
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
		for(int i = 0; i < gt.vc.size(); i++) vc.push_back(std::move(gt.vc[i]));
		//vc.insert(vc.end(), gt.vc.begin(), gt.vc.end());
		index[k].second = vc.size();
		//gt.vc.clear();
	}

	bridge_solver br(gx, vc, cfg, length_low, length_high);
	br.build_phase_set(px);

	// refine splice graph and phasing paths
	map<int32_t, int32_t> smap, tmap;
	group_start_boundaries(gx, smap, cfg.max_group_boundary_distance);
	group_end_boundaries(gx, tmap, cfg.max_group_boundary_distance);
	px.project_boundaries(smap, tmap);

	refine_splice_graph(gx);
	//keep_surviving_edges(gx, cfg.min_splicing_count);

	// construct hyper-set
	hyper_set hx(gx, px);
	hx.filter_nodes(gx);

	/*
	printf("---- parent\n");
	gx.print();
	hx.print_nodes();
	printf("====\n");
	*/

	// assemble combined graph
	gx.gid = cb.gid;
	cfg.algo = "single";
	scallop sx(gx, hx, &cfg);
	sx.assemble();

	for(int k = 0; k < sx.trsts.size(); k++)
	{
		transcript &t = sx.trsts[k];
		t.RPKM = 0;
		if(t.exons.size() <= 1) continue;
		index_transcript(mt, t);
		//t.write(cout);
		if(cfg.merge_intersection == false) vt.push_back(t);
	}

	// process each individual graph
	for(int i = 0; i < gv.size(); i++)
	{
		// process unbridged reads
		combined_graph g1;		
		for(int k = index[i].first; k < index[i].second; k++)
		{
			if(br.opt[k].type < 0) continue;
			g1.append(vc[k], br.opt[k]);
		}

		vector<combined_graph*> gv1;
		gv1.push_back(&g1);
		gv1.push_back(gv[i]);

		combined_graph cb1;
		cb1.combine(gv1);
		cb1.copy_meta_information(*(gv[i]));

		splice_graph gr;
		cb1.build_splice_graph(gr);
		gr.build_vertex_index();

		refine_splice_graph(gr);
		//keep_surviving_edges(gr, rs, cfg.min_splicing_count);

		gr.gid = gv[i]->gid;

		hyper_set hs(gr, cb1.ps);
		hs.filter_nodes(gr);

		/*
		printf("---- child %d\n", i);
		gr.print();
		hs.print_nodes();
		printf("====\n");
		*/

		cfg.algo = "single";
		scallop sc(gr, hs, &cfg);
		sc.assemble();

		for(int k = 0; k < sc.trsts.size(); k++)
		{
			transcript &t = sc.trsts[k];
			t.RPKM = 0;
			if(t.exons.size() <= 1) continue;
			//t.write(cout);
			bool b = query_transcript(mt, t);
			if(b == false) index_transcript(mt, t);
			if(b || cfg.merge_intersection == false) vt.push_back(t);
		}
	}

	printf("assemble combined-graph %s, %lu children, %lu assembled transcripts\n", cb.gid.c_str(), gv.size(), vt.size());
	return 0;
}

int assemble_cluster(vector<combined_graph*> gv, int instance, map< size_t, vector<transcript> > &trsts, mutex &mylock, const parameters &cfg)
{
	assert(gv.size() >= 2);
	int subindex = 0;
	vector<transcript> vt;
	for(int i = 0; i < gv.size(); i++)
	{
		for(int j = i + 1; j < gv.size(); j++)
		{
			vector<combined_graph*> gv1;
			gv1.push_back(gv[i]);
			gv1.push_back(gv[j]);
			assemble_cluster(gv1, instance, subindex++, vt, cfg);
		}
	}

	mylock.lock();
	for(int k = 0; k < vt.size(); k++)
	{
		index_transcript(trsts, vt[k]);
	}
	mylock.unlock();

	for(int i = 0; i > gv.size(); i++)
	{
		gv[i]->clear();
	}

	return 0;
}

int index_transcript(map< size_t, vector<transcript> > &mt, const transcript &t)
{
	//if(t.exons.size() <= 1) t.write(fout);
	if(t.exons.size() <= 1) return 0;

	size_t h = t.get_intron_chain_hashing();
	if(mt.find(h) == mt.end())
	{
		vector<transcript> v;
		v.push_back(t);
		mt.insert(pair<size_t, vector<transcript> >(h, v));
	}
	else
	{
		mt[h].push_back(t);
	}
	return 0;
}

bool query_transcript(const map< size_t, vector<transcript> > &mt, const transcript &t)
{
	size_t h = t.get_intron_chain_hashing();
	map< size_t, vector<transcript> >::const_iterator it = mt.find(h);
	if(it == mt.end()) return false;

	const vector<transcript> &v = it->second;
	for(int k = 0; k < v.size(); k++)
	{
		if(v[k].strand != t.strand) continue;
		bool b = v[k].intron_chain_match(t);
		if(b == true) return true;
	}

	return false;
}
