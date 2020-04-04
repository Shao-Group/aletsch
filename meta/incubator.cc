#include "incubator.h"
#include "generator.h"
#include "filter.h"
#include "cluster.h"
#include "meta_config.h"
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

incubator::incubator(const config &c)
{
	cfg = c;
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
	ifstream fin(input_bam_list.c_str());
	if(fin.fail())
	{
		printf("cannot open input-bam-list-file %s\n", input_bam_list.c_str());
		exit(0);
	}

	mutex mylock;								// lock for trsts
	boost::asio::thread_pool pool(max_threads); // thread pool

	char line[102400];
	while(fin.getline(line, 10240, '\n'))
	{
		string file(line);
		if(file.size() == 0) continue;
		boost::asio::post(pool, [this, &mylock, file]{ generate_single(file, this->groups, mylock, this->g2g, this->cfg); });
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
	boost::asio::thread_pool pool(max_threads); // thread pool

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
	boost::asio::thread_pool pool(max_threads); // thread pool
	mutex mylock;								// lock for trsts

	int instance = 0;
	for(int i = 0; i < groups.size(); i++)
	{
		for(int k = 0; k < groups[i].clusters.size(); k++)
		{
			vector<combined_graph> &gset = groups[i].gset;
			set<int> &cluster = groups[i].clusters[k];
			boost::asio::post(pool, [&gset, &cluster, instance, &mylock]{ assemble_single(gset, cluster, instance, this->trsts, mylock, this->cfg); });
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

	time_t mytime = time(NULL);
	printf("finish filtering and reporting final assembled transcripts, %s\n", ctime(&mytime));

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

int generate_single(const string &file, vector<combined_group> &gv, mutex &mylock, vector< map<string, int> > &m, const config &cfg)
{	
	vector<combined_graph> v;
	config c = cfg;
	c.input_file = file;
	generator gt(v, c);
	gt.resolve();

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

	time_t mytime = time(NULL);
	printf("finish processing individual sample %s, %s", file.c_str(), ctime(&mytime));
	return 0;
}

int assemble_single(const vector<combined_graph> &gset, const set<int> &cluster, int instance, map< size_t, vector<transcript> > &trsts, mutex &mylock, const config &cfg1)
{
	config cfg = cfg1;

	// setting up names
	char name[10240];
	sprintf(name, "instance.%d.0", instance);
	cb.gid = name;

	for(int j = 0; j < cb.children.size(); j++)
	{
		sprintf(name, "instance.%d.%d", instance, j + 1);
		cb.children[j].gid = name;
	}

	int z = 0;
	map< size_t, vector<transcript> > mt;
	vector<transcript> vt;

	set<int32_t> ps = cb.get_reliable_splices(min_supporting_samples, 99999);

	// rebuild splice graph and phasing paths
	splice_graph gx;
	phase_set px;
	cb.build_splice_graph(gx);
	cb.build_phase_set(px);

	// collect and bridge all unbridged pairs
	vector<pereads_cluster> preads;
	vector<PI> index(cb.children.size());
	int length_low = 999;
	int length_high = 0;
	for(int k = 0; k < cb.children.size(); k++)
	{
		combined_graph &gt = cb.children[k];
		if(gt.sp.insertsize_low < length_low) length_low = gt.sp.insertsize_low;
		if(gt.sp.insertsize_high > length_high) length_high = gt.sp.insertsize_high;
		index[k].first = preads.size();
		preads.insert(preads.end(), gt.preads.begin(), gt.preads.end());
		index[k].second = preads.size();
		gt.preads.clear();
	}

	bridge_solver br(gx, preads);
	br.length_low = length_low;
	br.length_high = length_high;
	br.build_phase_set(px);

	// refine splice graph and phasing paths
	map<int32_t, int32_t> smap, tmap;
	group_start_boundaries(gx, smap, max_group_boundary_distance);
	group_end_boundaries(gx, tmap, max_group_boundary_distance);
	px.project_boundaries(smap, tmap);

	refine_splice_graph(gx);
	keep_surviving_edges(gx, min_splicing_count);

	// construct hyper-set
	hyper_set hx(gx, px);

	// assemble parent graph
	gx.gid = cb.gid;
	cfg.algo = "single";
	scallop sx(gx, hx, &cfg);
	sx.assemble();

	for(int k = 0; k < sx.trsts.size(); k++)
	{
		transcript &t = sx.trsts[k];
		t.RPKM = 0;
		if(t.exons.size() <= 1) continue;
		z++;
		index_transcript(mt, t);
		if(merge_intersection == false || cb.children.size() == 0)
		{
			vt.push_back(t);
		}
	}

	for(int i = 0; i < cb.children.size(); i++)
	{
		// process unbridged reads
		phase_set pps;
		for(int k = index[i].first; k < index[i].second; k++)
		{
			if(br.opt[k].type < 0) continue;
			add_phase_from_bridge_path(preads[k], br.opt[k], pps);
			// TODO add junction to this child with br.opt[k]
			build_phase_from_bridge_path(gx, preads[k], br.opt[k], pps);
		}
		cb.children[i].combine_extra_phase_set(pps);

		/*
		if(exon_chains.size() >= 1)
		{
			mylock.lock();
			for(int k = 0; k < weights.size(); k++)
			{
				printf("extra read %d: weight = %d, list = ( ", k, weights[k]);
				printv(exon_chains[k]);
				printf(")\n");
			}
			cb.children[i].print(i);
			printf("-----\n");
			cb.children[i].combine_extra_bridged_reads(exon_chains, weights);
			cb.children[i].print(i);
			printf("=====\n");
			mylock.unlock();
		}
		*/

		splice_graph gr;
		phase_set ps;
		cb.children[i].resolve(gr, hs, ub);

		hyper_set hs;
		/*
		printf("-----\n");
		gr.print();
		hs.print_nodes();
		*/

		refine_splice_graph(gr);
		keep_surviving_edges(gr, ps, min_splicing_count);
		hs.filter_nodes(gr);

		gr.gid = cb.children[i].gid;

		cfg.algo = "single";
		scallop sc(gr, hs, &cfg);
		sc.assemble();

		int z1 = 0;
		int z2 = 0;
		for(int k = 0; k < sc.trsts.size(); k++)
		{
			transcript &t = sc.trsts[k];
			t.RPKM = 0;
			if(t.exons.size() <= 1) continue;

			z1++;
			bool b = query_transcript(mt, t);
			if(b == false) index_transcript(mt, t);
			if(b == true) z2++;
			if(b || merge_intersection == false)
			{
				vt.push_back(t);
			}
		}
	}

	printf("assemble combined-graph %s, %lu children, %lu assembled transcripts\n", cb.gid.c_str(), cb.children.size(), vt.size());

	if(vt.size() >= 1)
	{
		mylock.lock();
		for(int k = 0; k < vt.size(); k++)
		{
			index_transcript(trsts, vt[k]);
		}
		mylock.unlock();
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
