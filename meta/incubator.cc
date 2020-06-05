/*
Part of meta-scallop
(c) 2020 by Mingfu Shao, The Pennsylvania State University
See LICENSE for licensing.
*/

#include "incubator.h"
#include "generator.h"
#include "assembler.h"
#include "filter.h"
#include "cluster.h"
#include "parameters.h"
#include "scallop.h"
#include "graph_reviser.h"
#include "bridge_solver.h"
#include "essential.h"
#include "constants.h"
#include "previewer.h"

#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <thread>
#include <ctime>
#include <boost/asio/post.hpp>
#include <boost/asio/thread_pool.hpp>
#include <boost/pending/disjoint_sets.hpp>

incubator::incubator(const vector<parameters> &v)
	: params(v)
{
	meta_gtf.open(params[DEFAULT].output_gtf_file.c_str());
	if(meta_gtf.fail())
	{
		printf("cannot open output-get-file %s\n", params[DEFAULT].output_gtf_file.c_str());
		exit(0);
	}
}


incubator::~incubator()
{
	meta_gtf.close();
}

int incubator::resolve()
{
	read_bam_list();
	init_samples();
	build_sample_index();

	time_t mytime;
	for(auto &x: sindex)
	{
		string chrm = x.first;
		tmerge.chrm = chrm;

		// test
		//if(chrm != "chr1") continue;

		mytime = time(NULL);
		printf("start processing chrm %s\n", chrm.c_str());

		mytime = time(NULL);
		printf("step 1: generate graphs for individual bam/sam files, %s", ctime(&mytime));
		generate(x.second);

		mytime = time(NULL);
		printf("step 2: merge splice graphs, %s", ctime(&mytime));
		merge();

		mytime = time(NULL);
		printf("step 3: assemble merged splice graphs, %s", ctime(&mytime));
		assemble();

		groups.clear();

		mytime = time(NULL);
		printf("step 4: rearrange transcript sets, %s", ctime(&mytime));
		rearrange();

		mytime = time(NULL);
		printf("step 5: postprocess and write assembled transcripts, %s", ctime(&mytime));
		postprocess();

		mytime = time(NULL);
		printf("finish processing chrm %s\n", chrm.c_str());
	}

	close_samples();
	return 0;
}

int incubator::read_bam_list()
{
	ifstream fin(params[DEFAULT].input_bam_list.c_str());
	if(fin.fail())
	{
		printf("cannot open input-bam-list-file %s\n", params[DEFAULT].input_bam_list.c_str());
		exit(0);
	}

	char line[102400];
	while(fin.getline(line, 102400, '\n'))
	{
		if(strlen(line) <= 0) continue;
		stringstream sstr(line);
		sample_profile sp;
		char align_file[10240];
		char index_file[10240];
		char type[10240];
		sstr >> align_file >> index_file >> type;
		sp.align_file = align_file;
		sp.index_file = index_file;
		if(string(type) == "paired_end") sp.data_type = PAIRED_END;
		if(string(type) == "single_end") sp.data_type = SINGLE_END;
		if(string(type) == "pacbio_ccs") sp.data_type = PACBIO_CCS;
		if(string(type) == "pacbio_sub") sp.data_type = PACBIO_SUB;
		if(string(type) == "ont") sp.data_type = ONT;
		assert(sp.data_type != DEFAULT);
		sp.sample_id = samples.size();
		samples.push_back(sp);
	}
	return 0;
}

int incubator::init_samples()
{
	boost::asio::thread_pool pool(params[DEFAULT].max_threads);
	for(int i = 0; i < samples.size(); i++)
	{
		sample_profile &sp = samples[i];
		boost::asio::post(pool, [this, &sp]{ this->init_sample(sp); });
	}
	pool.join();
	return 0;
}

int incubator::build_sample_index()
{
	sindex.clear();
	for(int i = 0; i < samples.size(); i++)
	{
		sample_profile &sp = samples[i];
		for(int k = 0; k < sp.hdr->n_targets; k++)
		{
			string chrm(sp.hdr->target_name[k]);
			if(sindex.find(chrm) == sindex.end())
			{
				vector<PI> v;
				v.push_back(PI(i, k));
				sindex.insert(make_pair(chrm, v));
			}
			else
			{
				sindex[chrm].push_back(PI(i, k));
			}
		}
	}
	return 0;
}

int incubator::init_sample(sample_profile &sp)
{
	previewer pre(params[sp.data_type], sp);
	pre.infer_library_type();
	if(sp.data_type == PAIRED_END) pre.infer_insertsize();

	if(params[DEFAULT].output_bridged_bam_dir != "") 
	{
		sp.open_bridged_bam(params[DEFAULT].output_bridged_bam_dir);
	}

	if(params[DEFAULT].output_gtf_dir != "") 
	{
		sp.open_individual_gtf(params[DEFAULT].output_gtf_dir);
	}

	sp.open_align_file();
	sp.build_index_iterators();

	return 0;
}

int incubator::close_samples()
{
	for(int i = 0; i < samples.size(); i++)
	{
		if(params[DEFAULT].output_bridged_bam_dir != "") samples[i].close_bridged_bam();
		if(params[DEFAULT].output_gtf_dir != "") samples[i].close_individual_gtf();
		samples[i].destroy_index_iterators();
		samples[i].close_align_file();
	}
	return 0;
}

int incubator::generate(const vector<PI> &v)
{
	if(v.size() == 0) return 0;

	int num_threads = params[DEFAULT].max_threads;
	if(params[DEFAULT].single_sample_multiple_threading == true) num_threads = params[DEFAULT].max_threads / 2;
	
	boost::asio::thread_pool pool(num_threads);			// thread pool
	mutex mylock;										// lock for 

	for(int i = 0; i < v.size(); i++)
	{
		int sid = v[i].first;
		int tid = v[i].second;
		sample_profile &sp = samples[sid];
		boost::asio::post(pool, [this, &mylock, &sp, tid]{ this->generate(sp, tid, mylock); });
	}

	pool.join();
	//print_groups();

	return 0;
}

int incubator::merge()
{
	boost::asio::thread_pool pool(params[DEFAULT].max_threads); // thread pool

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
	boost::asio::thread_pool pool(params[DEFAULT].max_threads);
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

	return 0;
}

int incubator::rearrange()
{
	tmerge.mt.clear();
	for(int i = 0; i < tsets.size(); i++)
	{
		transcript_set &t = tsets[i];
		assert(t.chrm == tmerge.chrm);
		tmerge.add(t, 2, TRANSCRIPT_COUNT_ADD_COVERAGE_ADD);
	}
	tsets.clear();
	return 0;
}

int incubator::postprocess()
{
	vector<transcript> v = tmerge.get_transcripts(2);

	// warning: ts contains mixed strands
	//cluster cs(v, cfg);
	//cs.solve();

	filter ft(v, /*cs.cct,*/ params[DEFAULT]);
	//ft.join_single_exon_transcripts();
	ft.filter_length_coverage();

	stringstream ss;
	vector<vector<int>> vv(samples.size());
	for(int i = 0; i < ft.trs.size(); i++)
	{
		transcript &t = ft.trs[i];
		t.write(ss);
		pair<bool, trans_item> p = tmerge.query(t);
		if(p.first == false) continue;
		for(auto &k : p.second.samples)
		{
			if(k < 0) continue;
			assert(k >= 0 && k < vv.size());
			vv[k].push_back(i);
		}
	}

	const string &s = ss.str();
	meta_gtf.write(s.c_str(), s.size());

	boost::asio::thread_pool pool(params[DEFAULT].max_threads);
	for(int i = 0; i < vv.size(); i++)
	{
		const vector<transcript> &z = ft.trs;
		const vector<int> &v = vv[i];
		boost::asio::post(pool, [this, i, &z, &v]{ this->write_individual_gtf(i, z, v); });
	}
	pool.join();

	return 0;
}

int incubator::generate(sample_profile &sp, int tid, mutex &mylock)
{	
	vector<combined_graph> v;
	vector<transcript> trsts;
	generator gt(sp, v, trsts, params[sp.data_type], tid);
	gt.resolve();

	assert(trsts.size() == 0);
	//save_transcripts(trsts, sp.sample_id, mylock);

	mylock.lock();
	for(int k = 0; k < v.size(); k++)
	{
		bool found = false;
		for(int i = 0; i < groups.size(); i++)
		{
			if(groups[i].chrm != v[k].chrm) continue;
			if(groups[i].strand != v[k].strand) continue;
			groups[i].add_graph(std::move(v[k]));
			found = true;
			break;
		}

		if(found == false)
		{
			combined_group gp(v[k].chrm, v[k].strand, params[DEFAULT]);
			gp.add_graph(std::move(v[k]));
			groups.push_back(std::move(gp));
		}
	}
	mylock.unlock();
	printf("finish processing tid = %d of sample %s\n", tid, sp.align_file.c_str());
	return 0;
}

int incubator::assemble(vector<combined_graph*> gv, int instance, mutex &mylock)
{
	if(gv.size() == 0) return 0;

	transcript_set ts;

	assembler asmb(params[DEFAULT]);
	asmb.assemble(gv, 0, instance, ts, samples);

	move_transcript_set(ts, mylock);
	for(int i = 0; i > gv.size(); i++) gv[i]->clear();

	return 0;
}

int incubator::write_individual_gtf(int id, const vector<transcript> &vt, const vector<int> &v)
{
	assert(id >= 0 && id < samples.size());

	stringstream ss;
	for(int i = 0; i < v.size(); i++)
	{
		int k = v[i];
		const transcript &t = vt[k];
		t.write(ss);
	}

	const string &s = ss.str();

	sample_profile &sp = samples[id];
	sp.gtf_lock.lock();
	sp.individual_gtf->write(s.c_str(), s.size());
	sp.gtf_lock.unlock();

	return 0;
}

int incubator::move_transcript_set(transcript_set &ts, mutex &mylock)
{
	if(ts.mt.size() == 0) return 0;
	mylock.lock();
	tsets.push_back(std::move(ts));
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
