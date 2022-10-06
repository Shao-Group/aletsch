/*
Part of aletsch
(c) 2020 by Mingfu Shao, The Pennsylvania State University
See LICENSE for licensing.
*/

#include "incubator.h"
#include "generator.h"
#include "assembler.h"
#include "filter.h"
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

incubator::incubator(vector<parameters> &v)
	: params(v), tpool(params[DEFAULT].max_threads),
	tmerge("", params[DEFAULT].min_single_exon_clustering_overlap) 
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

	init_bundle_groups();

	// get max region

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

		int max_region = 0;
		for(auto &z: x.second)
		{
			if(max_region < samples[z.first].start1[z.second].size()) 
				max_region = samples[z.first].start1[z.second].size();
		}

		for(int k = 0; k < max_region; k++)
		{
			mytime = time(NULL);
			printf("processing chrm %s, region %d, max-region = %d, %s", chrm.c_str(), k, max_region, ctime(&mytime));

			//boost::asio::post(tpool, [this, chrm, k]{ 
			this->generate_merge_assemble(chrm, k);
			//});
		}

		break;

		// TODO
		/*
		mytime = time(NULL);
		printf("postprocess and write assembled transcripts for chrm %s, %s", chrm.c_str(), ctime(&mytime));
		rearrange();
		postprocess();

		mytime = time(NULL);
		printf("finish processing chrm %s, %s\n", chrm.c_str(), ctime(&mytime));
		*/
	}

	tpool.join();

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

int incubator::init_bundle_groups()
{
	for(auto &x: sindex)
	{
		string chrm = x.first;
		int max_region = 0;
		for(auto &z: x.second)
		{
			if(max_region < samples[z.first].start1[z.second].size()) 
				max_region = samples[z.first].start1[z.second].size();
		}

		for(int k = 0; k < max_region; k++)
		{
			bundle_group g1(chrm, '+', k, params[DEFAULT], tpool);
			bundle_group g2(chrm, '-', k, params[DEFAULT], tpool);
			bundle_group g3(chrm, '.', k, params[DEFAULT], tpool);
			grps.push_back(bundle_group(chrm, '+', k, params[DEFAULT], tpool));
			grps.push_back(bundle_group(chrm, '-', k, params[DEFAULT], tpool));
			grps.push_back(bundle_group(chrm, '.', k, params[DEFAULT], tpool));
		}
	}
	return 0;
}

int incubator::get_bundle_group(string chrm, int rid)
{
	for(int k = 0; k < grps.size(); k++)
	{
		if(grps[k].chrm != chrm) continue;
		if(grps[k].rid != rid) continue;
		return k;
	}
	return -1;
}

int incubator::generate_merge_assemble(string chrm, int rid)
{
	if(sindex.find(chrm) == sindex.end()) return 0;
	const vector<PI> &v = sindex[chrm];
	if(v.size() == 0) return 0;

	mutex group_lock;
	vector<mutex> sample_locks(v.size());
	for(int k = 0; k < sample_locks.size(); k++) sample_locks[k].lock();

	for(int i = 0; i < v.size(); i++)
	{
		int sid = v[i].first;
		int tid = v[i].second;
		sample_profile &sp = samples[sid];
		mutex &sample_lock = sample_locks[i];

		boost::asio::post(tpool, [this, &group_lock, &sample_lock, &sp, chrm, tid, rid]{ 
			this->generate(sp, tid, rid, chrm, group_lock, sample_lock); 
		});
	}

	for(int k = 0; k < sample_locks.size(); k++) sample_locks[k].lock();

	int bi = this->get_bundle_group(chrm, rid);
	for(int i = 0; i < 3; i++)
	{
		bundle_group &g = this->grps[bi + i];
		g.resolve(); 
		this->assemble(g, rid, i);
	}

	//for(int k = 0; k < sample_locks.size(); k++) sample_locks[k].unlock();

	return 0;
}

int incubator::generate(sample_profile &sp, int tid, int rid, string chrm, mutex &group_lock, mutex &sample_lock)
{	
	if(rid >= sp.start1.size()) return 0;

	vector<bundle> v;
	transcript_set ts(chrm, params[DEFAULT].min_single_exon_clustering_overlap);
	generator gt(sp, v, ts, params[sp.data_type], tid, rid);
	gt.resolve();
	save_transcript_set(ts, tlock);

	group_lock.lock();
	int bi = get_bundle_group(chrm, rid);
	assert(bi != -1);
	for(int k = 0; k < v.size(); k++)
	{
		if(v[k].strand == '+') grps[bi + 0].gset.push_back(std::move(v[k]));
		if(v[k].strand == '-') grps[bi + 1].gset.push_back(std::move(v[k]));
		if(v[k].strand == '.') grps[bi + 2].gset.push_back(std::move(v[k]));
	}
	group_lock.unlock();

	printf("finish processing tid = %d, rid = %d, of sample %s\n", tid, rid, sp.align_file.c_str());

	sample_lock.unlock();
	return 0;
}

int incubator::assemble(bundle_group &g, int rid, int gid)
{
	int instance = 0;
	vector<bool> vb(g.gset.size(), false);
	for(int k = 0; k < g.gvv.size(); k++)
	{
		const vector<int> &v = g.gvv[k];
		if(v.size() == 0) continue;
		vector<bundle*> gv;
		for(int j = 0; j < v.size(); j++)
		{
			gv.push_back(&(g.gset[v[j]]));
			assert(vb[v[j]] == false);
			vb[v[j]] = true;
		}
		boost::asio::post(tpool, [this, gv, instance, rid, gid]{ 
				//this->assemble(gv, instance, mylock, pool); 
				//transcript_set ts(gv.front()->chrm, params[DEFAULT].min_single_exon_clustering_overlap);
				assembler asmb(params[DEFAULT], this->tspool, this->tlock, this->tpool, rid, gid, instance);
				asmb.resolve(gv);
		});
		instance++;
	}
	return 0;
}

int incubator::rearrange()
{
	// filtering with count
	/*
	boost::asio::thread_pool pool(params[DEFAULT].max_threads);
	for(int i = 0; i < tspool.tsets.size(); i++)
	{
		transcript_set &t = tspool.tsets[i];
		assert(t.chrm == tmerge.chrm);
		boost::asio::post(pool, [&t]{ t.filter(2); });
	}
	pool.join();
	*/

	// random sort
	std::random_shuffle(tspool.tsets.begin(), tspool.tsets.end());

	// merge
	mutex mylock;
	boost::asio::thread_pool pool2(params[DEFAULT].max_threads);

	int t = params[DEFAULT].max_threads;
	if(t <= 0) t = 1;
	int n = ceil(1.0 * tspool.tsets.size() / t);

	tmerge.mt.clear();
	for(int i = 0; i < t; i++)
	{
		int a = (i + 0) * n;
		int b = (i + 1) * n;
		if(b >= tspool.tsets.size()) b = tspool.tsets.size();
		boost::asio::post(pool2, [this, &mylock, a, b]{ 
				transcript_set ts(this->tmerge.chrm, params[DEFAULT].min_single_exon_clustering_overlap);
				for(int k = a; k < b; k++) ts.add(this->tspool.tsets[k], TRANSCRIPT_COUNT_ADD_COVERAGE_ADD);
				mylock.lock();
				this->tmerge.add(ts, TRANSCRIPT_COUNT_ADD_COVERAGE_ADD);
				mylock.unlock();
			});
	}
	pool2.join();

	tspool.tsets.clear();
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
			//if(v[k].count <= 1) continue;

			transcript &t = v[k].trst;

			if(verify_length_coverage(t, params[DEFAULT]) == false) continue;
			if(verify_exon_length(t, params[DEFAULT]) == false) continue;

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

int incubator::write_individual_gtf(int id, const vector<transcript> &vt, const vector<int> &ct, const vector<pair<int, double>> &v)
{
	assert(id >= 0 && id < samples.size());

	stringstream ss;
	for(int i = 0; i < v.size(); i++)
	{
		int k = v[i].first;
		const transcript &t = vt[k];
		double cov2 = v[i].second;

		if(t.exons.size() == 1 && cov2 < params[DEFAULT].min_single_exon_individual_coverage) continue;

		t.write(ss, cov2, ct[k]);
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
	tspool.tsets.push_back(ts);
	mylock.unlock();
	return 0;
}

int incubator::print_groups(const vector<bundle_group> &groups)
{
	for(int k = 0; k < groups.size(); k++)
	{
		printf("group %d (chrm = %s, strand = %c) contains %lu graphs (%lu merged graphs)\n", k, groups[k].chrm.c_str(), groups[k].strand, groups[k].gset.size(), groups[k].gvv.size());
		//groups[k].print();
	}
	return 0;
}
