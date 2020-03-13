#include "incubator.h"
#include "generator.h"
#include "filter.h"
#include "cluster.h"
#include "meta_config.h"
#include "scallop.h"
#include "hyper_graph.h"
#include "graph_revise.h"
#include "boost/pending/disjoint_sets.hpp"
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <thread>
#include <boost/asio/thread_pool.hpp>
#include <boost/asio/post.hpp>

incubator::incubator(const config &c)
{
	cfg = c;
	g2g.resize(3);
}

int incubator::resolve()
{
	printf("loading bam/sam files ...\n");
	load();

	printf("merge splice graphs ...\n");
	merge();

	printf("assemble merged splice graphs ...\n");
	assemble();

	printf("filtering and output gtf files ...\n");
	postprocess();

	return 0;
}

int incubator::load()
{
	ifstream fin(input_bam_list.c_str());
	if(fin.fail())
	{
		printf("cannot open input-bam-list-file %s\n", input_bam_list.c_str());
		exit(0);
	}

	vector< vector<string> > files(2);

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
		threads.emplace_back(load_multiple, files[k], std::ref(groups), std::ref(mylock), std::ref(g2g), std::ref(cfg));
	}
	for(int k = 0; k < threads.size(); k++)
	{
		threads[k].join();
	}

	print_groups();
	return 0;
}

int incubator::merge()
{
	boost::asio::thread_pool pool(max_threads); // thread pool

	for(int k = 0; k < groups.size(); k++)
	{
		combined_group &gp = groups[k];
		boost::asio::post(pool, [&gp]{ gp.resolve(); });
	}

	pool.join();
	print_groups();

	return 0;	
}

int incubator::assemble()
{
	boost::asio::thread_pool pool(max_threads); // thread pool
	mutex mylock;								// lock for trsts

	int instance = 0;
	for(int i = 0; i < groups.size(); i++)
	{
		for(int k = 0; k < groups[i].mset.size(); k++)
		{
			combined_graph &cb = groups[i].mset[k];
			boost::asio::post(pool, [&cb, this, instance, &mylock]{ assemble_single(cb, instance, this->trsts, mylock, this->cfg); });
			instance++;
		}
	}
	pool.join();
	return 0;
}

int incubator::postprocess()
{
	ofstream fout(output_gtf_file.c_str());
	if(fout.fail())
	{
		printf("cannot open output-get-file %s\n", output_gtf_file.c_str());
		exit(0);
	}

	boost::asio::thread_pool pool(max_threads); // thread pool
	mutex mylock;								// lock for trsts

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
						cluster cs(v, &(this->cfg));
						cs.solve();

						filter ft(cs.cct, &(this->cfg));
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

	boost::asio::thread_pool pool(max_threads); // thread pool
	mutex mylock;								// lock for trsts
	for(int k = 0; k < groups.size(); k++)
	{
		combined_group &gp = groups[k];
		int os = offset[k];
		boost::asio::post(pool, [&gp, &fout, &mylock, os]{ gp.write(mylock, fout, os); });
	}
	pool.join();

	fout.flush();
	fout.close();
	return 0;
}

int incubator::print_groups()
{
	for(int k = 0; k < groups.size(); k++)
	{
		printf("group %d (chrm = %s, strand = %c) contains %lu graphs (%lu merged graphs)\n", k, groups[k].chrm.c_str(), groups[k].strand, groups[k].gset.size(), groups[k].mset.size());
		/*
		for(int j = 0; j < groups[k].gset.size(); j++)
		{
			groups[k].gset[j].print(j);
		}
		*/
	}
	return 0;
}

int load_multiple(const vector<string> &files, vector<combined_group> &gv, mutex &mylock, vector< map<string, int> > &m, const config &cfg)
{	
	vector<combined_graph> v;
	for(int k = 0; k < files.size(); k++) 
	{
		//printf("load file %s\n", files[k].c_str());
		//load_single(files[k], v);
		config c = cfg;
		c.input_file = files[k];
		generator gt(v, c);
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

int assemble_single(combined_graph &cb, int instance, map< size_t, vector<transcript> > &trsts, mutex &mylock, const config &cfg1)
{
	config cfg = cfg1;

	char name[10240];
	sprintf(name, "instance.%d.0", instance);
	merged_graph cm;
	cm.build(cb, name);
	if(merge_intersection == true) cm.parent = false;
	else cm.parent = true;

	vector<merged_graph> children(cb.children.size());
	for(int j = 0; j < cb.children.size(); j++)
	{
		sprintf(name, "instance.%d.%d", instance, j + 1);
		children[j].build(cb.children[j], name);
		children[j].parent = false;
	}

	//if(cm.num_combined <= 0) return 0;

	int z = 0;
	map< size_t, vector<transcript> > mt;
	vector<transcript> vt;

	set<int32_t> ps = cm.get_reliable_splices(min_supporting_samples, 99999);
	set<int32_t> sb = cm.get_reliable_start_boundaries(min_supporting_samples, 99999);
	set<int32_t> tb = cm.get_reliable_end_boundaries(min_supporting_samples, 99999);
	set<int32_t> aj = cm.get_reliable_adjacencies(min_supporting_samples, min_splicing_count);
	set<PI32> rs = cm.get_reliable_junctions(min_supporting_samples, min_splicing_count);

	cm.solve();

	keep_surviving_edges(cm.gr, min_splicing_count);
	//keep_surviving_edges(cm.gr, ps, min_splicing_count);

	cm.hs.filter_nodes(cm.gr);

	//cm.hs.print_nodes();
	//cm.hs.filter();

	//cm.gr.print();
	//cm.hs.print_nodes();

	string gid = cm.gr.gid;
	cm.gr.gid = gid + ".0";

	cfg.algo = "single";
	scallop sm(cm.gr, cm.hs, &cfg);

	//sm.gr.print();
	//sm.hs.print_nodes();

	sm.assemble();

	//for(int k = 0; k < sm.paths.size(); k++) sm.paths[k].print(k);


	for(int k = 0; k < sm.trsts.size(); k++)
	{
		transcript &t = sm.trsts[k];
		t.RPKM = 0;
		if(t.exons.size() <= 1) continue;
		z++;
		index_transcript(mt, t);
		//t.write(cout);
		//if(merge_intersection == false || children.size() == 0) index_transcript(trsts, t);
		if(merge_intersection == false || children.size() == 0)
		{
			vt.push_back(t);
			//mylock.lock();
			//index_transcript(trsts, t);
			//mylock.unlock();
		}
	}

	//printf("--------\n");

	for(int i = 0; i < children.size(); i++)
	{
		merged_graph &cb = children[i];
		cb.solve();

		keep_surviving_edges(cb.gr, ps, min_splicing_count);
		//filter_graph(cb.gr, ps, aj, sb, tb, min_splicing_count);
		//keep_surviving(cb.gr, ps, aj, sb, tb, min_splicing_count);
		//filter_junctions(cb.gr, ps, min_splicing_count);
		//filter_start_boundaries(cb.gr, sb, min_splicing_count);
		//filter_end_boundaries(cb.gr, tb, min_splicing_count);

		cb.hs.filter_nodes(cb.gr);

		//cb.hs.print_nodes();

		//cb.build_phasing_paths(cm.paths);
		//cb.hx.filter_nodes(cb.gr);
		//cb.hx.filter();

		/*
		printf("--------------\n");
		hyper_graph hg(cb.hs.nodes);
		hg.keep_maximal_nodes();
		hg.build_overlap_index();
		hg.print_nodes();
		hg.print_index();
		hg.align_paths(cb.hx.nodes);
		printf("--------------\n\n");
		*/

		//cb.hs.extend(cb.hx);
		//cb.hs.filter();

		cb.gr.gid = gid + "." + tostring(i + 1);

		cfg.algo = "single";
		//scallop sc(cb.gr, cb.hs, cb.hx);
		scallop sc(cb.gr, cb.hs, &cfg);
		sc.assemble();

		//for(int k = 0; k < sc.paths.size(); k++) sc.paths[k].print(k);


		int z1 = 0;
		int z2 = 0;
		for(int k = 0; k < sc.trsts.size(); k++)
		{
			transcript &t = sc.trsts[k];
			t.RPKM = 0;
			//if(t.exons.size() <= 1) t.write(fout);
			if(t.exons.size() <= 1) continue;

			//t.write(cout);

			z1++;
			bool b = query_transcript(mt, t);
			if(b == false) index_transcript(mt, t);
			if(b == true) z2++;
			//if(b == false) 
			if(b || merge_intersection == false)
			{
				vt.push_back(t);
				//mylock.lock();
				//index_transcript(trsts, t);
				//mylock.unlock();
			}
		}

		printf("combined-graph %s, %d transcripts, %d child %s, %d -> %d transcripts\n", cm.gid.c_str(), z, i, cb.gr.gid.c_str(), z1, z2);
	}

	if(vt.size() >= 1)
	{
		mylock.lock();
		for(int k = 0; k < vt.size(); k++)
		{
			index_transcript(trsts, vt[k]);
		}
		mylock.unlock();
	}

	//printf("=====\n");
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

	printf("hash fail");
	return false;
}
