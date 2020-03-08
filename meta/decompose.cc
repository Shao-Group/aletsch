#include "filter.h"
#include "cluster.h"
#include "config.h"
#include "scallop.h"
#include "merged_graph.h"
#include "hyper_graph.h"
#include "graph_revise.h"
#include "decompose.h"

#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <string>
#include <mutex>
#include <boost/asio/thread_pool.hpp>
#include <boost/asio/post.hpp>

int assemble()
{
	ifstream fin(graph_file.c_str());
	if(fin.fail())
	{
		printf("could not open file %s\n", graph_file.c_str());
		exit(0);
	}

	char line[10240];
	char gid[10240];
	char chrm[10240];
	char mark[1024];
	char strand[1024];
	int num_combined;

	ofstream fout(output_file.c_str());
	if(fout.fail()) return 0;

	map< size_t, vector<transcript> > trsts;
	boost::asio::thread_pool pool(max_threads); // thread pool
	mutex mylock;								// lock for trsts

	merged_graph mgraph;
	vector<merged_graph> children;

	int index = -1;
	while(fin.getline(line, 10240, '\n'))
	{
		if(line[0] != '#') continue;
		stringstream sstr(line);
		num_combined = 0;
		sstr >> mark >> gid >> chrm >> strand >> num_combined;

		//printf("num_combined = %d, index = %d\n", num_combined, index);

		if(index <= 0) 
		{
			if(index == 0) 
			{
				//assemble(mgraph, children, trsts, mylock);
				//boost::asio::post(pool, [index]{ test(index); });
				boost::asio::post(pool, [mgraph, children, &trsts, &mylock]{ assemble(mgraph, children, trsts, mylock); });
			}

			mgraph.clear();
			children.clear();

			mgraph.gid = gid;

			if(merge_intersection == true) mgraph.parent = false;
			else mgraph.parent = true;
			//mgraph.parent = true;

			mgraph.build(fin, gid, chrm, strand[0], num_combined);
			index = num_combined;
			if(index <= 1) index = 0;
		}
		else
		{
			merged_graph cb;
			cb.parent = false;
			cb.build(fin, gid, chrm, strand[0], num_combined);
			children.push_back(cb);
			index--;
		}
	}

	pool.join();

	boost::asio::thread_pool pool2(max_threads); // thread pool

	typedef map<size_t, vector<transcript> >::iterator MIT;
	index = 0;
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

		boost::asio::post(pool2, [m1, m2, &fout, &mylock]
				{ 
					stringstream ss;
					for(MIT x = m1; x != m2; x++)
					{
						vector<transcript> &v = x->second;
						cluster cs(v);
						cs.solve();

						filter ft(cs.cct);
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

	pool2.join();

	fin.close();
	fout.close();
	return 0;
}

int assemble(merged_graph cm, vector<merged_graph> children, map< size_t, vector<transcript> > &trsts, mutex &mylock)
{
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

	algo = "single";
	scallop sm(cm.gr, cm.hs);

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

		algo = "single";
		//scallop sc(cb.gr, cb.hs, cb.hx);
		scallop sc(cb.gr, cb.hs);
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
