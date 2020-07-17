/*
Part of aletsch
(c) 2020 by Mingfu Shao, The Pennsylvania State University
See LICENSE for licensing.
*/

#include "incubator.h"
#include "generator.h"
#include "assembler.h"
#include "filter.h"
#include "cluster.h"
#include "parameters.h"
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

incubator::incubator(vector<parameters> &v)
	: params(v), tmerge("", params[DEFAULT].min_single_exon_clustering_overlap)
{
	if(params[DEFAULT].profile_only == true) return;
	meta_gtf.open(params[DEFAULT].output_gtf_file.c_str());
	if(meta_gtf.fail())
	{
		printf("cannot open output-get-file %s\n", params[DEFAULT].output_gtf_file.c_str());
		exit(0);
	}
}

incubator::~incubator()
{
	if(params[DEFAULT].profile_only == true) return;
	meta_gtf.close();
}

int incubator::resolve()
{
	read_bam_list();
	init_samples();

	if(params[DEFAULT].profile_only == true) return 0;

	build_sample_index();

	time_t mytime;
	for(auto &x: sindex)
	{
		string chrm = x.first;
		tmerge.chrm = chrm;

		// test
		/*
		size_t found = chrm.find("_");
		if(found == string::npos) continue;
		if(chrm != "chr11_KI270721v1_random") continue;
		if(chrm == "chr1") continue;
		if(chrm == "chr10") continue;
		if(chrm == "chr11") continue;
		if(chrm == "chr12") continue;
		if(chrm == "chr13") continue;
		if(chrm == "chr14") continue;
		if(chrm == "chr11_KI270721v1_random") continue;
		if(chrm == "chr14_GL000009v2_random") continue;
		*/

		mytime = time(NULL);
		printf("start processing chrm %s, %s", chrm.c_str(), ctime(&mytime));

		mytime = time(NULL);
		printf("step 1: generate graphs for individual bam/sam files, %s", ctime(&mytime));
		generate(chrm);

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
		printf("finish processing chrm %s, %s\n", chrm.c_str(), ctime(&mytime));
	}

	free_samples();
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
		sample_profile sp(samples.size());
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
		assert(sp.sample_id == samples.size());
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
		boost::asio::post(pool, [this, &sp] 
		{
				const parameters &cfg = this->params[sp.data_type];

				if(cfg.profile_only == true)
				{
					previewer pre(cfg, sp);
					pre.infer_library_type();
					if(sp.data_type == PAIRED_END) pre.infer_insertsize();
					if(cfg.profile_dir != "") sp.save_profile(cfg.profile_dir);
					return;
				}

				if(cfg.profile_dir != "")
				{
					sp.load_profile(cfg.profile_dir);
				}
				else
				{
					previewer pre(cfg, sp);
					pre.infer_library_type();
					if(sp.data_type == PAIRED_END) pre.infer_insertsize();
				}

				sp.read_index_iterators(); 
				string bdir = cfg.output_bridged_bam_dir;
				if(bdir != "") sp.init_bridged_bam(bdir);
		});
	}
	pool.join();
	return 0;
}

int incubator::free_samples()
{
	for(int i = 0; i < samples.size(); i++) 
	{
		samples[i].free_index_iterators();
		//samples[i].free_align_headers();
	}
	return 0;
}

int incubator::build_sample_index()
{
	set<string> ss;
	if(params[DEFAULT].chrm_list_file != "")
	{
		ifstream fin(params[DEFAULT].chrm_list_file.c_str());
		if(fin.fail()) printf("cannot open chrm list file\n");
		if(fin.fail()) exit(0);
		char line[10240];
		while(fin.getline(line, 10240, '\n'))
		{
			if(string(line) == "") continue;
			ss.insert(string(line));
		}
		fin.close();
	}

	if(params[DEFAULT].chrm_list_string != "")
	{
		vector<string> v = split_string(params[DEFAULT].chrm_list_string, ",");
		for(int i = 0; i < v.size(); i++)
		{
			if(v[i] == "") continue;
			ss.insert(v[i]);
		}
	}

	sindex.clear();
	for(int i = 0; i < samples.size(); i++)
	{
		sample_profile &sp = samples[i];
		sp.open_align_file();
		for(int k = 0; k < sp.hdr->n_targets; k++)
		{
			string chrm(sp.hdr->target_name[k]);
			if(ss.size() >= 1 && ss.find(chrm) == ss.end()) continue;

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
		sp.close_align_file();
	}
	return 0;
}

int incubator::generate(string chrm)
{
	if(sindex.find(chrm) == sindex.end()) return 0;
	const vector<PI> &v = sindex[chrm];
	if(v.size() == 0) return 0;

	int num_threads = params[DEFAULT].max_threads;
	boost::asio::thread_pool pool(num_threads);			// thread pool
	mutex mylock;										// lock for 

	for(int i = 0; i < v.size(); i++)
	{
		int sid = v[i].first;
		int tid = v[i].second;
		sample_profile &sp = samples[sid];
		boost::asio::post(pool, [this, &mylock, &sp, chrm, tid]{ this->generate(sp, tid, chrm, mylock); });
	}

	pool.join();
	print_groups();

	return 0;
}

int incubator::merge()
{
	for(int k = 0; k < groups.size(); k++) groups[k].resolve();
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
		vector<bool> vb(groups[i].gset.size(), false);
		for(int k = 0; k < groups[i].gvv.size(); k++)
		{
			const vector<int> &v = groups[i].gvv[k];
			if(v.size() == 0) continue;
			vector<combined_graph*> gv;
			for(int j = 0; j < v.size(); j++)
			{
				gv.push_back(&(groups[i].gset[v[j]]));
				assert(vb[v[j]] == false);
				vb[v[j]] = true;
			}
			boost::asio::post(pool, [this, gv, instance, &mylock]{ this->assemble(gv, instance, mylock); });
			instance++;
		}
	}
	pool.join();

	return 0;
}

int incubator::rearrange()
{
	// filtering with count
	boost::asio::thread_pool pool(params[DEFAULT].max_threads);
	for(int i = 0; i < tsets.size(); i++)
	{
		transcript_set &t = tsets[i];
		assert(t.chrm == tmerge.chrm);
		boost::asio::post(pool, [&t]{ t.filter(2); });
	}
	pool.join();

	// random sort
	std::random_shuffle(tsets.begin(), tsets.end());

	// merge
	mutex mylock;
	boost::asio::thread_pool pool2(params[DEFAULT].max_threads);

	int t = params[DEFAULT].max_threads;
	if(t <= 0) t = 1;
	int n = ceil(1.0 * tsets.size() / t);

	tmerge.mt.clear();
	for(int i = 0; i < t; i++)
	{
		int a = (i + 0) * n;
		int b = (i + 1) * n;
		if(b >= tsets.size()) b = tsets.size();
		boost::asio::post(pool2, [this, &mylock, a, b]{ 
				transcript_set ts(this->tmerge.chrm, params[DEFAULT].min_single_exon_clustering_overlap);
				for(int k = a; k < b; k++) ts.add(this->tsets[k], TRANSCRIPT_COUNT_ADD_COVERAGE_ADD);
				mylock.lock();
				this->tmerge.add(ts, TRANSCRIPT_COUNT_ADD_COVERAGE_ADD);
				mylock.unlock();
			});
	}
	pool2.join();

	tsets.clear();
	return 0;
}

int incubator::postprocess()
{
	stringstream ss;
	vector<transcript> vt;
	vector<int> ct;
	vector<vector<pair<int, double>>> vv(samples.size());
	for(auto &it : tmerge.mt)
	{
		auto &v = it.second;
		for(int k = 0; k < v.size(); k++)
		{
			if(v[k].count <= 1) continue;

			transcript &t = v[k].trst;
			if(verify_length_coverage(t, params[DEFAULT]) == false) continue;

			t.write(ss, -1, v[k].samples.size());
			vt.push_back(t);
			ct.push_back(v[k].samples.size());
			
			for(auto &p : v[k].samples)
			{
				int j = p.first;
				double w = p.second;
				if(j < 0 || j >= vv.size()) continue;
				vv[j].push_back(make_pair(vt.size() - 1, w));
			}
		}
	}

	const string &s = ss.str();
	meta_gtf.write(s.c_str(), s.size());

	if(params[DEFAULT].output_gtf_dir != "")
	{
		boost::asio::thread_pool pool(params[DEFAULT].max_threads);
		for(int i = 0; i < vv.size(); i++)
		{
			const vector<int> &c = ct;
			const vector<transcript> &z = vt;
			const vector<pair<int, double>> &v = vv[i];
			boost::asio::post(pool, [this, i, &z, &c, &v]{ this->write_individual_gtf(i, z, c, v); });
		}
		pool.join();
	}

	return 0;
}

int incubator::generate(sample_profile &sp, int tid, string chrm, mutex &mylock)
{	
	vector<combined_graph> v;
	transcript_set ts(chrm, params[DEFAULT].min_single_exon_clustering_overlap);
	generator gt(sp, v, ts, params[sp.data_type], tid);
	gt.resolve();
	save_transcript_set(ts, mylock);

	mylock.lock();
	for(int k = 0; k < v.size(); k++)
	{
		bool found = false;
		for(int i = 0; i < groups.size(); i++)
		{
			if(groups[i].chrm != v[k].chrm) continue;
			if(groups[i].strand != v[k].strand) continue;
			groups[i].gset.push_back(std::move(v[k]));
			found = true;
			break;
		}

		if(found == false)
		{
			graph_group gp(v[k].chrm, v[k].strand, params[DEFAULT]);
			gp.gset.push_back(std::move(v[k]));
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

	transcript_set ts(gv.front()->chrm, params[DEFAULT].min_single_exon_clustering_overlap);

	//printf("assemble instance %d with %lu graphs\n", instance, gv.size());
	//for(int k = 0; k < gv.size(); k++) gv[k]->print(k);

	assembler asmb(params[DEFAULT]);
	asmb.assemble(gv, 0, instance, ts, samples);

	save_transcript_set(ts, mylock);
	for(int i = 0; i > gv.size(); i++) gv[i]->clear();

	return 0;
}

int incubator::write_individual_gtf(int id, const vector<transcript> &vt, const vector<int> &ct, const vector<pair<int, double>> &v)
{
	assert(id >= 0 && id < samples.size());

	stringstream ss;
	for(int i = 0; i < v.size(); i++)
	{
		int k = v[i].first;
		transcript t = vt[k];
		double cov1 = v[i].second;

		if(t.exons.size() == 1 && t.coverage < params[DEFAULT].min_single_exon_individual_coverage) continue;

		t.write(ss, cov1, ct[k]);
	}

	const string &s = ss.str();

	sample_profile &sp = samples[id];
	sp.gtf_lock.lock();
	sp.open_individual_gtf(params[DEFAULT].output_gtf_dir);
	sp.individual_gtf->write(s.c_str(), s.size());
	sp.close_individual_gtf();
	sp.gtf_lock.unlock();

	return 0;
}

int incubator::save_transcript_set(const transcript_set &ts, mutex &mylock)
{
	if(ts.mt.size() == 0) return 0;
	mylock.lock();
	tsets.push_back(ts);
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
