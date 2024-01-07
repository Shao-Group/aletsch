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
	: params(v), tpool(params[DEFAULT].max_threads), group_size(params[DEFAULT].max_threads), gmutex(99999), tmutex(99999)
{
	if(params[DEFAULT].profile_only == true) return;
	meta_gtf.open(params[DEFAULT].output_gtf_file.c_str(), std::ofstream::out | std::ofstream::app);
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
	init_transcript_sets();

	//for(int k = 0; k < samples.size(); k++) samples[k].open_align_file();

	time_t mytime;
	for(auto &x: sindex)
	{
		string chrm = x.first;
		int m = ceil(get_max_region(chrm) * 1.0 / group_size);
		//int m = get_max_region(chrm);
		for(int k = 0; k < m; k++)
		{
			generate_merge_assemble(chrm, k);
		}
	}
	tpool.join();

	mytime = time(NULL);
	printf("postprocess and write assembled transcripts, %s", ctime(&mytime));
	postprocess();

	mytime = time(NULL);
	printf("free samples, %s", ctime(&mytime));
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
		sample_profile sp(samples.size(), params[DEFAULT].region_partition_length);
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
	if(samples.size() <= 0) return 0;

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
				//if(cfg.profile_dir != "") sp.save_profile(cfg.profile_dir);
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
		});
	}
	pool.join();

	// post process previewed profile, to allow
	// for samples with fewer reads to borrow 
	// profiles from samples with more reads

	const parameters &cfg = this->params[DEFAULT];

	// profiles are given
	if(cfg.profile_only == false && cfg.profile_dir != "") return 0;

	int b = 0;
	for(int i = 1; i < samples.size(); i++)
	{
		if(samples[i].spn > samples[b].spn)
		{
			b = i;
		}
		else if(samples[i].spn == samples[b].spn && samples[i].insert_total > samples[b].insert_total) 
		{
			b = i;
		}
	}

	for(int i = 0; i < samples.size(); i++)
	{
		// copy profile from 
		if(samples[i].insert_total < cfg.min_preview_spliced_reads)
		{
			samples[i].insertsize_ave = samples[b].insertsize_ave; 
			samples[i].insertsize_std = samples[b].insertsize_std; 
			samples[i].insertsize_low = samples[b].insertsize_low; 
			samples[i].insertsize_high = samples[b].insertsize_high; 
		}

		if(samples[i].spn < cfg.min_preview_spliced_reads)
		{
			samples[i].library_type = samples[b].library_type; 
			samples[i].bam_with_xs = samples[b].bam_with_xs; 
		}
				
		if(cfg.profile_dir != "") samples[i].save_profile(cfg.profile_dir);
	}

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

int incubator::get_chrm_index(string chrm, int sid)
{
	assert(sindex.find(chrm) != sindex.end());

	for(auto &z: sindex[chrm])
	{
		if(z.first == sid) return z.second;
	}
	return -1;
}

int incubator::get_max_region(string chrm)
{
	assert(sindex.find(chrm) != sindex.end());

	int max_region = 0;
	for(auto &z: sindex[chrm])
	{
		//printf("get max-region: chrm = %s, z.first = %d, z.second = %d, samples[z.first].start1[z.second].size = %lu\n", chrm.c_str(), z.first, z.second, samples[z.first].start1[z.second].size());
		if(max_region < samples[z.first].start1[z.second].size()) 
			max_region = samples[z.first].start1[z.second].size();
	}
	return max_region;
}

int incubator::init_bundle_groups()
{
	grps.clear();
	for(auto &z : sindex)
	{
		string chrm = z.first;
		//int m = ceil(get_max_region(chrm) * 1.0 / group_size);
		int m = get_max_region(chrm);
		for(int k = 0; k < m; k++)
		{
			grps.push_back(bundle_group(chrm, '+', k, params[DEFAULT], sindex));
			grps.push_back(bundle_group(chrm, '-', k, params[DEFAULT], sindex));
			grps.push_back(bundle_group(chrm, '.', k, params[DEFAULT], sindex));
		}
	}
	return 0;
}

int incubator::init_transcript_sets()
{
	tts.clear();
	for(auto &z : sindex)
	{
		string chrm = z.first;
		transcript_set t1(chrm, -9, params[DEFAULT].min_single_exon_clustering_overlap);
		transcript_set t2(chrm, -9, params[DEFAULT].min_single_exon_clustering_overlap);
		transcript_set t3(chrm, -9, params[DEFAULT].min_single_exon_clustering_overlap);
		tts.insert(make_pair(make_pair(chrm, '+'), t1));
		tts.insert(make_pair(make_pair(chrm, '-'), t2));
		tts.insert(make_pair(make_pair(chrm, '.'), t3));
	}
	return 0;
}

int incubator::get_bundle_group(string chrm, int rid)
{
	for(int k = 0; k < grps.size(); k++)
	{
		//if(grps[k].chrm == chrm && grps[k].gid * group_size <= rid && rid < (grps[k].gid + 1) * group_size) return k;
		if(grps[k].chrm == chrm && grps[k].rid == rid) return k;
	}
	return -1;
}

int incubator::generate_merge_assemble(string chrm, int gid)
{
	if(sindex.find(chrm) == sindex.end()) return 0;
	const vector<PI> &v = sindex[chrm];
	if(v.size() == 0) return 0;

	vector<mutex> curlocks(v.size() * group_size);
	for(int k = 0; k < curlocks.size(); k++) curlocks[k].lock();

	for(int j = 0; j < group_size; j++)
	{
		for(int i = 0; i < v.size(); i++)
		{
			int sid = v[i].first;
			int tid = v[i].second;
			//sample_profile &sp = samples[sid];

			int rid = gid * group_size + j;
			mutex &curlock = curlocks[i * group_size + j];

			time_t mytime = time(NULL);
			//printf("generate chrm %s, gid = %d, rid = %d, %s", chrm.c_str(), gid, rid, ctime(&mytime));
			boost::asio::post(this->tpool, [this, &curlock, sid, chrm, tid, rid]{ 
					this->generate(sid, tid, rid, chrm, curlock); 
			});
		}
	}

	//for(int k = 0; k < curlocks.size(); k++) curlocks[k].lock();

	// print start/end positions
	/*
	for(int i = 0; i < v.size(); i++)
	{
		int sid = v[i].first;
		int tid = v[i].second;
		for(int j = 0; j < group_size; j++)
		{
			int rid = gid * group_size + j;
			if(rid >= samples[sid].start1[tid].size()) continue;
			if(rid >= samples[sid].start2[tid].size()) continue;
			printf("sample %d, tid = %d, rid = %d, strand +, expected end = %d, actual end = %d\n", sid, tid, rid, samples[sid].start1[tid][rid] + samples[sid].region_partition_length, samples[sid].end1[tid][rid]);
			printf("sample %d, tid = %d, rid = %d, strand -, expected end = %d, actual end = %d\n", sid, tid, rid, samples[sid].start2[tid][rid] + samples[sid].region_partition_length, samples[sid].end2[tid][rid]);
		}
	}
	*/

	vector<bool> posted(group_size, false);
	while(true)
	{
		for(int j = 0; j < group_size; j++)
		{
			if(posted[j] == true) continue;

			bool succeed0 = true;
			bool succeed1 = true;
			vector<bool> ck0(v.size(), false);
			vector<bool> ck1(v.size(), false);

			if(j >= 1)
			{
				for(int i = 0; i < v.size(); i++) 
				{
					ck0[i] = curlocks[i * group_size + j - 1].try_lock();
					if(ck0[i] == false) succeed0 = false;
				}
			}

			for(int i = 0; i < v.size(); i++) 
			{
				ck1[i] = curlocks[i * group_size + j].try_lock();
				if(ck1[i] == false) succeed1 = false;
			}

			if(succeed0 && succeed1)
			{
				int rid = gid * group_size + j;
				int bi = this->get_bundle_group(chrm, rid);
				if(bi >= 0)
				{
					for(int i = 0; i < 3; i++)
					{
						bundle_group &g = this->grps[bi + i];
						time_t mytime = time(NULL);
						//printf("assemble chrm %s, gid = %d, rid = %d, bi = %d, %s", chrm.c_str(), gid, rid, bi, ctime(&mytime));
						boost::asio::post(this->tpool, [this, &g, bi, rid, i]{ 
								g.resolve(); 
								this->assemble(g, rid, i);
								g.clear();
						});
					}
				}
				posted[j] = true;
			}

			for(int i = 0; i < v.size(); i++) 
			{
				if(j >= 1 && ck0[i] == true) curlocks[i * group_size + j - 1].unlock();
				if(j >= 0 && ck1[i] == true) curlocks[i * group_size + j - 0].unlock();
			}
		}

		bool all_posted = true;
		for(int j = 0; j < posted.size(); j++)
		{
			if(posted[j] == false) all_posted = false;
		}
		
		if(all_posted == true) break;
	}

	//for(int k = 0; k < samples.size(); k++) samples[k].close_align_file();
	return 0;
}

int incubator::generate(int sid, int tid, int rid, string chrm, mutex &curlock)
{	
	sample_profile &sp = samples[sid];
	int cid = get_chrm_index(chrm, sid);

	//printf("sp.start1[cid].size = %lu, sid = %d, rid = %d, tid = %d, cid = %d, chrm = %s\n", sp.start1[cid].size(), sid, rid, tid, cid, chrm.c_str());

	if(rid >= sp.start1[cid].size()) 
	{
		//printf("unlock rid = %d\n", rid);
		curlock.unlock();
		return 0;
	}

	vector<bundle> v;
	//transcript_set ts(chrm, params[DEFAULT].min_single_exon_clustering_overlap);

	//printf("in generating sid = %d, tid = %d, rid = %d, chrm = %s, cid = %d, regions = %lu\n", sid, tid, rid, chrm.c_str(), cid, sp.start1[cid].size());
	int bi = get_bundle_group(chrm, rid);
	assert(bi != -1);

	generator gt(sp, v, params[sp.data_type], tid, rid);
	gt.resolve();

	//save_transcript_set(ts, tlock);

	gmutex[bi + 0].lock();
	for(int k = 0; k < v.size(); k++)
	{
		if(v[k].strand == '+' && v[k].splices.size() >= 1) grps[bi + 0].gset.emplace_back(std::move(v[k]));
	}
	gmutex[bi + 0].unlock();

	gmutex[bi + 1].lock();
	for(int k = 0; k < v.size(); k++)
	{
		if(v[k].strand == '-' && v[k].splices.size() >= 1) grps[bi + 1].gset.emplace_back(std::move(v[k]));
	}
	gmutex[bi + 1].unlock();

	gmutex[bi + 2].lock();
	for(int k = 0; k < v.size(); k++)
	{
		if(v[k].strand == '.' && v[k].splices.size() >= 1) grps[bi + 2].gset.emplace_back(std::move(v[k]));
	}
	gmutex[bi + 2].unlock();

	curlock.unlock();

	mutex mtx;
	transcript_set ts0(chrm, rid, params[DEFAULT].min_single_exon_clustering_overlap);
	transcript_set ts1(chrm, rid, params[DEFAULT].min_single_exon_clustering_overlap);
	transcript_set ts2(chrm, rid, params[DEFAULT].min_single_exon_clustering_overlap);

	int index = 0;
	int cnt0 = 0, cnt1 = 0, cnt2 = 0;
	for(int k = 0; k < v.size(); k++)
	{
		if(v[k].splices.size() >= 1) continue;

		assert(v[k].splices.size() <= 0);

		if(v[k].strand == '+')
		{
			assembler asmb(params[DEFAULT], ts0, mtx, rid, sid, index++);
			asmb.assemble(v[k]);
			cnt0++;
		}
		if(v[k].strand == '-')
		{
			assembler asmb(params[DEFAULT], ts1, mtx, rid, sid, index++);
			asmb.assemble(v[k]);
			cnt1++;
		}
		if(v[k].strand == '.')
		{
			assembler asmb(params[DEFAULT], ts2, mtx, rid, sid, index++);
			asmb.assemble(v[k]);
			cnt2++;
		}
	}

	if(cnt0 >= 1)
	{
		//vector<transcript> v = ts0.get_transcripts(1);
		//for(int i = 0; i < v.size(); i++) v[i].write(cout);
		tmutex[bi + 0].lock();
		grps[bi + 0].num_assembled += cnt0;
		grps[bi + 0].tmerge.add(ts0, TRANSCRIPT_COUNT_ADD_COVERAGE_ADD);
		tmutex[bi + 0].unlock();
	}

	if(cnt1 >= 1)
	{
		tmutex[bi + 1].lock();
		grps[bi + 1].num_assembled += cnt1;
		grps[bi + 1].tmerge.add(ts1, TRANSCRIPT_COUNT_ADD_COVERAGE_ADD);
		tmutex[bi + 1].unlock();
	}

	if(cnt2 >= 1)
	{
		tmutex[bi + 2].lock();
		grps[bi + 2].num_assembled += cnt2;
		grps[bi + 2].tmerge.add(ts2, TRANSCRIPT_COUNT_ADD_COVERAGE_ADD);
		tmutex[bi + 2].unlock();
	}

	printf("finish generating tid = %d, rid = %d, of sample %s\n", tid, rid, sp.align_file.c_str());
	return 0;
}

int incubator::assemble(bundle_group &g, int rid, int gi)
{
	int instance = g.num_assembled + 1;
	vector<bool> vb(g.gset.size(), false);
	int sid = samples.size();
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
		assert(g.rid == rid);
		int bi = get_bundle_group(g.chrm, rid);
		mutex &mtx = tmutex[bi + gi];
		boost::asio::post(this->tpool, [this, &g, &mtx, gv, rid, sid, instance]{ 
				assembler asmb(params[DEFAULT], g.tmerge, mtx, rid, sid, instance);
				asmb.resolve(gv);
		});
		instance++;
	}
	return 0;
}

/*
int incubator::rearrange(transcript_set_pool &tsp, transcritp_set &tmerge)
{
	// random sort
	//std::random_shuffle(tspool.tsets.begin(), tspool.tsets.end());

	// merge
	mutex mylock;

	int t = params[DEFAULT].max_threads;
	if(t <= 0) t = 1;
	int n = ceil(1.0 * tspool.tsets.size() / t);

	tmerge.mt.clear();
	boost::asio::thread_pool pool(params[DEFAULT].max_threads);
	for(int i = 0; i < t; i++)
	{
		int a = (i + 0) * n;
		int b = (i + 1) * n;
		if(b >= tspool.tsets.size()) b = tspool.tsets.size();
		boost::asio::post(pool, [this, &mylock, a, b]{ 
				transcript_set ts(this->tmerge.chrm, params[DEFAULT].min_single_exon_clustering_overlap);
				for(int k = a; k < b; k++) ts.add(this->tspool.tsets[k], TRANSCRIPT_COUNT_ADD_COVERAGE_ADD);
				mylock.lock();
				this->tmerge.add(ts, TRANSCRIPT_COUNT_ADD_COVERAGE_ADD);
				mylock.unlock();
		});
	}
	pool.join();
	tspool.clear();
	return 0;
}
*/

int incubator::postprocess()
{
	boost::asio::thread_pool pool1(params[DEFAULT].max_threads);
	for(auto &z: tts)
	{
		string chrm = z.first.first;
		char strand = z.first.second;
		transcript_set &ts = z.second;
		boost::asio::post(pool1, [this, chrm, strand, &ts] 
		{
			for(int k = 0; k < this->grps.size(); k++)
			{
				if(this->grps[k].chrm != chrm) continue;
				if(this->grps[k].strand != strand) continue;
				//printf("arrange chrm %s, strand %c, grp %d\n", chrm.c_str(), strand, k);
				ts.add(this->grps[k].tmerge, TRANSCRIPT_COUNT_ADD_COVERAGE_ADD);
				//this->grps[k].tmerge.clear();
			}
		});
	}
	pool1.join();

	// extend samples
	sample_profile sn(samples.size(), samples[0].region_partition_length);
	samples.emplace_back(std::move(sn));

	// write gtf
	boost::asio::thread_pool pool2(params[DEFAULT].max_threads);
	boost::asio::post(pool2, [this] { this->write_combined_gtf(); });
	for(int i = 0; i < samples.size(); i++)
	{
		boost::asio::post(pool2, [this, i]{ this->write_individual_gtf(i); });
	}
	pool2.join();
	return 0;
}

int incubator::write_combined_gtf()
{
	for(auto &z: tts)
	{
		string chrm = z.first.first;
		char strand = z.first.second;
		const transcript_set &tm = z.second;

		stringstream ss;
		for(auto &it : tm.mt)
		{
			auto &v = it.second;
			for(int k = 0; k < v.size(); k++)
			{
				const transcript &t = v[k].trst;

				//if(verify_length_coverage(t, params[DEFAULT]) == false) continue;
				//if(verify_exon_length(t, params[DEFAULT]) == false) continue;
				assert(v[k].samples.size() == t.count2);
				t.write(ss, -1, v[k].samples.size());

				//if(t.exons.size() > 1) t.write_features(-1);
				//Only output novel transcripts in merged graph
				
				// FIXME
				if(t.exons.size() > 1 && t.count2 == 1 && v[k].samples.find(-1) != v[k].samples.end()) t.write_features(-1);
			}
		}
		const string &s = ss.str();
		meta_gtf.write(s.c_str(), s.size());
	}
	return 0;
}

int incubator::write_individual_gtf(int sid)
{
	sample_profile &sp = samples[sid];
	sp.gtf_lock.lock();
	sp.open_individual_gtf(params[DEFAULT].output_gtf_dir);

	//sp.individual_gtf->write(s.c_str(), s.size());

	for(auto &z: tts)
	{
		string chrm = z.first.first;
		char strand = z.first.second;
		const transcript_set &tm = z.second;

		stringstream ss;
		for(auto &it : tm.mt)
		{
			auto &v = it.second;
			for(int k = 0; k < v.size(); k++)
			{
				for(auto &p : v[k].samples)
				{
					int j = p.first;
					if(j == -1) j = samples.size() - 1;
					if(j != sid) continue;

					const transcript &t = p.second;

					assert(p.second.count2 == t.count2);
					assert(abs(p.second.coverage - t.coverage)<SMIN);

					if(t.exons.size() == 1 && t.cov2 < params[DEFAULT].min_single_exon_individual_coverage) continue;
					//t.write(*(sp.individual_gtf), t.cov2, t.count2);
					t.write(ss, t.cov2, t.count2);

					// FIXME
					//if(t.exons.size() > 1) t.write_features(sid);
				}
			}
		}
		const string &s = ss.str();
		sp.individual_gtf->write(s.c_str(), s.size());
	}

	sp.close_individual_gtf();
	sp.gtf_lock.unlock();

	return 0;
}

int incubator::write_individual_gtf(int id, const vector<transcript> &v)
{
	assert(id >= 0 && id < samples.size());

	stringstream ss;
	for(int i = 0; i < v.size(); i++)
	{
		const transcript &t = v[i];

		double cov2 = v[i].cov2;
		if(t.exons.size() == 1 && cov2 < params[DEFAULT].min_single_exon_individual_coverage) continue;

		int ct = v[i].count2;
		t.write(ss, cov2, ct);

        if(t.exons.size() > 1) t.write_features(id);
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

/*
int incubator::save_transcript_set(const transcript_set &ts, mutex &mylock)
{
	if(ts.mt.size() == 0) return 0;
	mylock.lock();
	tspool.tsets.push_back(ts);
	mylock.unlock();
	return 0;
}
*/

int incubator::print_groups(const vector<bundle_group> &groups)
{
	for(int k = 0; k < groups.size(); k++)
	{
		printf("group %d (chrm = %s, strand = %c) contains %lu graphs (%lu merged graphs)\n", k, groups[k].chrm.c_str(), groups[k].strand, groups[k].gset.size(), groups[k].gvv.size());
		//groups[k].print();
	}
	return 0;
}
