/*
Part of Coral
(c) 2019 by Mingfu Shao, The Pennsylvania State University.
See LICENSE for licensing.
*/

#include <cassert>
#include <cstdio>
#include <map>
#include <iomanip>
#include <fstream>

#include "bundle2.h"
#include "region.h"
#include "config.h"
#include "util.h"
#include "undirected_graph.h"
#include "bridger.h"

bundle2::bundle2(const bundle_base &bb)
	: bundle_base(bb)
{
	//compute_strand();
}

bundle2::~bundle2()
{}

int bundle2::build()
{
	build_junctions();
	extend_junctions();
	build_regions();

	//link_regions();

	align_hits_transcripts();
	index_references();

	build_fragments();
	group_fragments();
	bridge_fragments();

	return 0;
}

int bundle2::build_junctions()
{
	int min_max_boundary_quality = min_mapping_quality;
	map< int64_t, vector<int> > m;
	for(int i = 0; i < hits.size(); i++)
	{
		vector<int64_t> v = hits[i].spos;
		if(v.size() == 0) continue;

		//hits[i].print();
		for(int k = 0; k < v.size(); k++)
		{
			int64_t p = v[k];

			// DEBUG
			/*
			int32_t x1 = low32(p);
			int32_t x2 = high32(p);
			if(fabs(x1 - 9364768) <= 2 || fabs(x2 - 9364768) <=2)
			{
				printf("HIT ");
				hits[i].print();
			}
			*/
			if(m.find(p) == m.end())
			{
				vector<int> hv;
				hv.push_back(i);
				m.insert(pair< int64_t, vector<int> >(p, hv));
			}
			else
			{
				m[p].push_back(i);
			}
		}
	}

	map< int64_t, vector<int> >::iterator it;
	for(it = m.begin(); it != m.end(); it++)
	{
		vector<int> &v = it->second;
		if(v.size() < min_splice_boundary_hits) continue;

		int32_t p1 = high32(it->first);
		int32_t p2 = low32(it->first);

		int s0 = 0;
		int s1 = 0;
		int s2 = 0;
		int nm = 0;
		for(int k = 0; k < v.size(); k++)
		{
			hit &h = hits[v[k]];
			nm += h.nm;
			if(h.xs == '.') s0++;
			if(h.xs == '+') s1++;
			if(h.xs == '-') s2++;
		}

		//printf("junction: %s:%d-%d (%d, %d, %d) %d\n", chrm.c_str(), p1, p2, s0, s1, s2, s1 < s2 ? s1 : s2);

		junction jc(it->first, v.size());
		jc.nm = nm;
		if(s1 == 0 && s2 == 0) jc.strand = '.';
		else if(s1 >= 1 && s2 >= 1) jc.strand = '.';
		else if(s1 > s2) jc.strand = '+';
		else jc.strand = '-';
		junctions.push_back(jc);

		/*
		uint32_t max_qual = 0;
		for(int k = 0; k < v.size(); k++)
		{
			hit &h = hits[v[k]];
			if(h.qual > max_qual) max_qual = h.qual;
		}
		assert(max_qual >= min_max_boundary_quality);
		*/
	}
	return 0;
}

int bundle2::extend_junctions()
{
	map< int64_t, vector<int> > m;
	for(int i = 0; i < ref_trsts.size(); i++)
	{
		vector<PI32> v = ref_trsts[i].get_intron_chain();
		for(int k = 0; k < v.size(); k++)
		{
			assert(v[k].first < v[k].second);
			// TODO, TODO
			if(v[k].first <= lpos) continue;
			if(v[k].second >= rpos) continue;
			int64_t p = pack(v[k].first, v[k].second);

			if(m.find(p) == m.end())
			{
				vector<int> hv;
				hv.push_back(i);
				m.insert(pair< int64_t, vector<int> >(p, hv));
			}
			else
			{
				m[p].push_back(i);
			}
		}
	}

	map< int64_t, vector<int> >::iterator it;
	for(it = m.begin(); it != m.end(); it++)
	{
		vector<int> &v = it->second;

		int s0 = 0;
		int s1 = 0;
		int s2 = 0;
		for(int k = 0; k < v.size(); k++)
		{
			char c = ref_trsts[v[k]].strand;
			if(c == '.') s0++;
			if(c == '+') s1++;
			if(c == '-') s2++;
		}

		junction jc(it->first, 0 - v.size());
		jc.nm = 0;
		if(s1 == 0 && s2 == 0) jc.strand = '.';
		else if(s1 >= 1 && s2 >= 1) jc.strand = '.';
		else if(s1 > s2) jc.strand = '+';
		else jc.strand = '-';
		junctions.push_back(jc);
	}
	return 0;
}


int bundle2::build_regions()
{
	MPI s;
	s.insert(PI(lpos, START_BOUNDARY));
	s.insert(PI(rpos, END_BOUNDARY));
	for(int i = 0; i < junctions.size(); i++)
	{
		junction &jc = junctions[i];

		int32_t l = jc.lpos;
		int32_t r = jc.rpos;

		if(s.find(l) == s.end()) s.insert(PI(l, LEFT_SPLICE));
		else if(s[l] == RIGHT_SPLICE) s[l] = LEFT_RIGHT_SPLICE;

		if(s.find(r) == s.end()) s.insert(PI(r, RIGHT_SPLICE));
		else if(s[r] == LEFT_SPLICE) s[r] = LEFT_RIGHT_SPLICE;
	}

	vector<PPI> v(s.begin(), s.end());
	sort(v.begin(), v.end());

	regions.clear();
	for(int k = 0; k < v.size() - 1; k++)
	{
		int32_t l = v[k].first;
		int32_t r = v[k + 1].first;
		int ltype = v[k].second; 
		int rtype = v[k + 1].second; 

		if(ltype == LEFT_RIGHT_SPLICE) ltype = RIGHT_SPLICE;
		if(rtype == LEFT_RIGHT_SPLICE) rtype = LEFT_SPLICE;

		regions.push_back(region(l, r, ltype, rtype));
	}

	return 0;
}

int bundle2::link_regions()
{
	if(regions.size() == 0) return 0;

	MPI lm;
	MPI rm;
	for(int i = 0; i < regions.size(); i++)
	{
		int32_t l = regions[i].lpos;
		int32_t r = regions[i].rpos;

		//printf("region %d [%d, %d), lpos = %d, rpos = %d\n", i, l, r, lpos, rpos);

		assert(lm.find(l) == lm.end());
		assert(rm.find(r) == rm.end());
		lm.insert(PPI(l, i));
		rm.insert(PPI(r, i));
	}

	for(int i = 0; i < junctions.size(); i++)
	{
		junction &b = junctions[i];

		//printf("look for junction %d [%d, %d)\n", i, b.lpos, b.rpos);

		MPI::iterator li = rm.find(b.lpos);
		MPI::iterator ri = lm.find(b.rpos);

		assert(li != rm.end());
		assert(ri != lm.end());

		b.lregion = li->second;
		b.rregion = ri->second;
	}
	return 0;
}

int bundle2::align_hits_transcripts()
{
	map<int32_t, int> m;
	for(int k = 0; k < regions.size(); k++)
	{
		if(k >= 1) assert(regions[k - 1].rpos == regions[k].lpos);
		m.insert(pair<int32_t, int>(regions[k].lpos, k));
	}

	for(int i = 0; i < hits.size(); i++)
	{
		align_hit(m, hits[i], hits[i].vlist);
		hits[i].vlist = encode_vlist(hits[i].vlist);
	}

	ref_phase.resize(ref_trsts.size());
	for(int i = 0; i < ref_trsts.size(); i++)
	{
		align_transcript(m, ref_trsts[i], ref_phase[i]);
	}

	return 0;
}

int bundle2::align_hit(const map<int32_t, int> &m, const hit &h, vector<int> &vv)
{
	vv.clear();
	vector<int64_t> v;
	h.get_aligned_intervals(v);
	if(v.size() == 0) return 0;

	vector<PI> sp;
	sp.resize(v.size());

	int32_t p1 = high32(v.front());
	int32_t p2 = low32(v.back());

	sp[0].first = locate_region(p1);
	for(int k = 1; k < v.size(); k++)
	{
		p1 = high32(v[k]);

		map<int32_t, int>::const_iterator it = m.find(p1);

		assert(it != m.end());
		sp[k].first = it->second;
	}

	sp[sp.size() - 1].second = locate_region(p2 - 1);
	for(int k = 0; k < v.size() - 1; k++)
	{
		p2 = low32(v[k]);
		map<int32_t, int>::const_iterator it = m.find(p2);
		assert(it != m.end());
		sp[k].second = it->second - 1; 
	}

	for(int k = 0; k < sp.size(); k++)
	{
		assert(sp[k].first <= sp[k].second);
		if(k > 0) assert(sp[k - 1].second < sp[k].first);
		for(int j = sp[k].first; j <= sp[k].second; j++) vv.push_back(j);
	}

	return 0;
}

int bundle2::align_transcript(const map<int32_t, int> &m, const transcript &t, vector<int> &vv)
{
	vv.clear();
	int k1 = -1;
	int k2 = -1;
	for(int k = 0; k < t.exons.size(); k++)
	{
		if(t.exons[k].second > lpos)
		{
			k1 = k;
			break;
		}
	}
	for(int k = t.exons.size() - 1; k >= 0; k--)
	{
		if(t.exons[k].first < rpos)
		{
			k2 = k;
			break;
		}
	}

	if(k1 > k2) return 0;
	if(k1 == -1 || k2 == -1) return 0;

	vector<PI> sp;
	sp.resize(k2 + 1);

	int32_t p1 = t.exons[k1].first > lpos ? t.exons[k1].first : lpos;
	int32_t p2 = t.exons[k2].second < rpos ? t.exons[k2].second : rpos;

	sp[k1].first = locate_region(p1);
	for(int k = k1 + 1; k <= k2; k++)
	{
		p1 = t.exons[k].first;
		map<int32_t, int>::const_iterator it = m.find(p1);
		assert(it != m.end());
		sp[k].first = it->second;
	}

	sp[k2].second = locate_region(p2 - 1);
	for(int k = k1; k < k2; k++)
	{
		p2 = t.exons[k].second;
		map<int32_t, int>::const_iterator it = m.find(p2);
		assert(it != m.end());
		sp[k].second = it->second - 1; 
	}

	for(int k = k1; k <= k2; k++)
	{
		assert(sp[k].first <= sp[k].second);
		if(k > k1) assert(sp[k - 1].second < sp[k].first);
		for(int j = sp[k].first; j <= sp[k].second; j++) vv.push_back(j);
	}

	return 0;
}

int bundle2::index_references()
{
	ref_index.clear();
	ref_index.resize(regions.size());
	for(int k = 0; k < ref_phase.size(); k++)
	{
		vector<int> &v = ref_phase[k];
		for(int j = 0; j < v.size(); j++)
		{
			int x = v[j];
			ref_index[x].push_back(PI(k, j));
		}
	}
	return 0;
}

int bundle2::locate_region(int32_t x)
{
	if(regions.size() == 0) return -1;

	int k1 = 0;
	int k2 = regions.size();
	while(k1 < k2)
	{
		int m = (k1 + k2) / 2;
		region &r = regions[m];
		if(x >= r.lpos && x < r.rpos) return m;
		else if(x < r.lpos) k2 = m;
		else k1 = m;
	}
	return -1;
}

int bundle2::build_fragments()
{
	// TODO: parameters
	int32_t max_misalignment1 = 20;
	int32_t max_misalignment2 = 10;

	fragments.clear();
	if(hits.size() == 0) return 0;

	int max_index = hits.size() + 1;
	if(max_index > 1000000) max_index = 1000000;

	vector< vector<int> > vv;
	vv.resize(max_index);

	// first build index
	for(int i = 0; i < hits.size(); i++)
	{
		hit &h = hits[i];
		if(h.isize >= 0) continue;
		if(h.vlist.size() == 0) continue;

		// do not use hi; as long as qname, pos and isize are identical
		int k = (h.qhash % max_index + h.pos % max_index + (0 - h.isize) % max_index) % max_index;
		/*
		SI si(h.qname, h.hi);
		MSI &m = vv[k];
		assert(m.find(si) == m.end());
		m.insert(PSI(si, i));
		*/
		vv[k].push_back(i);
	}

	for(int i = 0; i < hits.size(); i++)
	{
		hit &h = hits[i];
		if(h.paired == true) continue;
		if(h.isize <= 0) continue;
		if(h.vlist.size() == 0) continue;

		int k = (h.qhash % max_index + h.mpos % max_index + h.isize % max_index) % max_index;

		/*
		h.print();
		for(int j = 0; j < vv[k].size(); j++)
		{
			hit &z = hits[vv[k][j]];
			printf(" ");
			z.print();
		}
		*/

		int x = -1;
		for(int j = 0; j < vv[k].size(); j++)
		{
			hit &z = hits[vv[k][j]];
			//if(z.hi != h.hi) continue;
			if(z.paired == true) continue;
			if(z.pos != h.mpos) continue;
			if(z.isize + h.isize != 0) continue;
			if(z.qhash != h.qhash) continue;
			if(z.qname != h.qname) continue;
			x = vv[k][j];
			break;
		}

		/*
		SI si(h.qname, h.hi);
		MSI::iterator it = vv[k].find(si);
		if(it == vv[k].end()) continue;
		int x = it->second;
		*/

		//printf("HIT: i = %d, x = %d, hits[i].vlist = %lu | ", i, x, hits[i].vlist.size(), hits[i].qname.c_str()); hits[i].print();

		if(x == -1) continue;
		if(hits[x].vlist.size() == 0) continue;

		fragment fr(&hits[i], &hits[x]);
		fr.lpos = h.pos;
		fr.rpos = hits[x].rpos;

		vector<int> v1 = decode_vlist(hits[i].vlist);
		vector<int> v2 = decode_vlist(hits[x].vlist);
		fr.k1l = fr.h1->pos - regions[v1.front()].lpos;
		fr.k1r = regions[v1.back()].rpos - fr.h1->rpos;
		fr.k2l = fr.h2->pos - regions[v2.front()].lpos;
		fr.k2r = regions[v2.back()].rpos - fr.h2->rpos;

		fr.b1 = true;
		if(v1.size() <= 1) 
		{
			fr.b1 = false;
		}
		else if(v1.size() >= 2 && v1[v1.size() - 2] == v1.back() - 1)
		{
			if(fr.h1->rpos - regions[v1.back()].lpos > max_misalignment1 + fr.h1->nm) fr.b1 = false;
		}
		else if(v1.size() >= 2 && v1[v1.size() - 2] != v1.back() - 1)
		{
			if(fr.h1->rpos - regions[v1.back()].lpos > max_misalignment2 + fr.h1->nm) fr.b1 = false;
		}

		fr.b2 = true;
		if(v2.size() <= 1)
		{
			fr.b2 = false;
		}
		else if(v2.size() >= 2 || v2[1] == v2.front() + 1)
		{
			if(regions[v2.front()].rpos - fr.h2->pos > max_misalignment1 + fr.h2->nm) fr.b2 = false;
		}
		else if(v2.size() >= 2 || v2[1] != v2.front() + 1)
		{
			if(regions[v2.front()].rpos - fr.h2->pos > max_misalignment2 + fr.h2->nm) fr.b2 = false;
		}

		fragments.push_back(fr);

		hits[i].paired = true;
		hits[x].paired = true;
	}

	//printf("total hits = %lu, total fragments = %lu\n", hits.size(), fragments.size());
	return 0;
}

int bundle2::group_fragments()
{
	if(fragments.size() == 0) return 0;

	sort(fragments.begin(), fragments.end(), compare_fragment);

	vector<fragment> ff;

	fragment fx = fragments[0];
	assert(fx.h1->vlist.size() >= 1);
	assert(fx.h2->vlist.size() >= 1);
	for(int k = 1; k < fragments.size(); k++)
	{
		fragment &fr = fragments[k];
		assert(fr.h1->vlist.size() >= 1);
		assert(fr.h2->vlist.size() >= 1);

		if(fx.equal(fr) == true)
		{
			fx.append(fr);
		}
		else
		{
			ff.push_back(fx);
			fx = fr;
		}
	}
	ff.push_back(fx);
	fragments = ff;

	//printf("grouped fragments = %lu\n", fragments.size());
	return 0;
}

int bundle2::bridge_fragments()
{
	bridger bdg(this);
	bdg.bridge();

	return 0;
}

int32_t bundle2::compute_aligned_length(int32_t k1l, int32_t k2r, const vector<int>& v)
{
	if(v.size() == 0) return 0;
	int32_t flen = 0;
	for(int i = 0; i < v.size(); i++)
	{
		int k = v[i];
		flen += regions[k].rpos - regions[k].lpos;
	}
	return flen - k1l - k2r;
}

int bundle2::build_hyper_nodes()
{
	map<vector<int>, int> m;

	for(int k = 0; k < fragments.size(); k++)
	{
		fragment &fr = fragments[k];

		assert(fr.h1->paired == true);
		assert(fr.h2->paired == true);

		if(fr.h1->bridged == false) continue;
		if(fr.h2->bridged == false) continue;

		assert(fr.paths.size() == 1);

		vector<int> &v = fr.paths[0].v;
		
		if(m.find(v) == m.end()) m.insert(pair<vector<int>, int>(v, 1));
		else m[v]++;
	}

	for(int k = 0; k < hits.size(); k++)
	{
		hit &h = hits[k];
		if(h.bridged == true) continue;

		vector<int> &v = h.vlist;
		
		if(m.find(v) == m.end()) m.insert(pair<vector<int>, int>(v, 1));
		else m[v]++;
	}

	hpnodes.clear();
	for(map<vector<int>, int>::iterator it = m.begin(); it != m.end(); it++)
	{
		vector<int> v = decode_vlist(it->first);
		int c = it->second;
		hpnodes.push_back(hyper_node(v, v, c));
	}

	return 0;
}

int bundle2::print(int index)
{
	printf("Bundle %d: ", index);

	// statistic xs
	int n0 = 0, np = 0, nq = 0;
	for(int i = 0; i < hits.size(); i++)
	{
		if(hits[i].xs == '.') n0++;
		if(hits[i].xs == '+') np++;
		if(hits[i].xs == '-') nq++;
	}

	printf("tid = %d, #hits = %lu, #fragments = %lu, #ref-trsts = %lu, range = %s:%d-%d, orient = %c (%d, %d, %d)\n",
			tid, hits.size(), fragments.size(), ref_trsts.size(), chrm.c_str(), lpos, rpos, strand, n0, np, nq);

	// print ref-trsts
	//for(int k = 0; k < ref_trsts.size(); k++) ref_trsts[k].write(cout);

	// print fragments 
	//for(int i = 0; i < fragments.size(); i++) fragments[i].print(i);

	/*
	// print fclusters
	for(int i = 0; i < fclusters.size(); i++) fclusters[i].print(i);
	*/

	if(verbose <= 1) return 0;

	// print junctions 
	for(int i = 0; i < junctions.size(); i++)
	{
		junctions[i].print(chrm, i);
	}

	// print hits
	for(int i = 0; i < hits.size(); i++) hits[i].print();

	// print regions
	for(int i = 0; i < regions.size(); i++)
	{
		regions[i].print(i);
	}

	// print junctions 
	for(int i = 0; i < junctions.size(); i++)
	{
		junctions[i].print(chrm, i);
	}

	printf("\n");

	return 0;
}

int bundle2::build_bam1_t(bam1_t &b1t, const vector<int> &v, const hit &h, int32_t p1, int32_t p2, char c)
{
	b1t.core = h;
	b1t.core.pos = p1;
	b1t.core.bin = 0;
	b1t.core.n_cigar = 0;
	b1t.core.l_qseq = 0;
	b1t.core.isize = p2 - p1;
	b1t.core.flag = h.flag;
	//b1t.core.flag -= b1t.core.flag & (0x1);

	b1t.m_data = b1t.core.l_qname + 4 * (v.size() * 2 - 1) + 7 * 3;
	b1t.data = new uint8_t[b1t.m_data];

	// copy qname
	b1t.l_data = 0;
	assert(h.qname.size() + b1t.core.l_extranul + 1 == b1t.core.l_qname);
	memcpy(b1t.data, h.qname.c_str(), h.qname.size());
	b1t.l_data += h.qname.length();
	for(int i = 0; i <= b1t.core.l_extranul; i++)
	{
		b1t.data[b1t.l_data] = 0;
		b1t.l_data++;
	}
	assert(b1t.l_data == b1t.core.l_qname);

	int skip = 0;
	// CIGAR
	int32_t x1 = p1;
	int32_t x2 = p1;
	for(int i = 0; i < v.size(); i++)
	{
		region &r = regions[v[i]];
		if(x1 < x2 && x2 < r.lpos)
		{
			add_cigar_match(b1t, x1, x2);
			add_cigar_skip(b1t, x2, r.lpos);
			skip++;
			x1 = r.lpos;
			x2 = (r.rpos < p2) ? r.rpos : p2;
		}
		else
		{
			x2 = (r.rpos < p2) ? r.rpos : p2;
		}
	}

	if(x1 < x2) add_cigar_match(b1t, x1, x2);

	if(c !='.') bam_aux_append(&(b1t), "XS", 'A', 1, (uint8_t*)(&c));
	if(h.hi != -1) bam_aux_append(&(b1t), "HI", 'i', 4, (uint8_t*)(&h.hi));
	if(h.nh != -1) bam_aux_append(&(b1t), "NH", 'i', 4, (uint8_t*)(&h.nh));

	return skip;
}

int bundle2::add_cigar_match(bam1_t &b1t, int32_t p1, int32_t p2)
{
	assert(p1 < p2);
	uint32_t c1 = p2 - p1;
	c1 = c1 << 4;
	c1 += BAM_CMATCH;
	memcpy(b1t.data + b1t.l_data, &c1, 4);
	b1t.l_data += 4;
	b1t.core.n_cigar++;
	return 0;
}

int bundle2::add_cigar_skip(bam1_t &b1t, int32_t p1, int32_t p2)
{
	assert(p1 < p2);
	uint32_t c1 = p2 - p1;
	c1 = c1 << 4;
	c1 += BAM_CREF_SKIP;
	memcpy(b1t.data + b1t.l_data, &c1, 4);
	b1t.l_data += 4;
	b1t.core.n_cigar++;
	return 0;
}

int bundle2::write_bam(BGZF *fout, map<PI, bam1_t> &bmap, map<int, bam1_t> &umap)
{
	for(int k = 0; k < fragments.size(); k++)
	{
		fragment &fr = fragments[k];

		assert(fr.h1->paired == true);
		assert(fr.h2->paired == true);

		if(fr.h1->bridged == false) continue;
		if(fr.h2->bridged == false) continue;

		assert(fr.paths.size() == 1);
		assert(fr.h1->qname == fr.h2->qname);

		hit *ht = fr.h1;
		hit *hx = fr.h2;
		while(ht != NULL && hx != NULL)
		{
			assert(ht->paired == true);
			assert(hx->paired == true);
			assert(ht->bridged == true);
			assert(hx->bridged == true);

			//if(ht->nh >= 2) bset.insert(ht->qname);

			char c = ht->xs;
			if(c == '.' && hx->xs != '.') c = hx->xs;

			//printf("hit %s ht = %c, hx = %c, c = %c\n", ht->qname.c_str(), ht->xs, hx->xs, c);
		
			bam1_t b1t;
			build_bam1_t(b1t, decode_vlist(fr.paths[0].v), *ht, fr.lpos, fr.rpos, c);
			b1t.core.mpos = 0;

			if(library_type == UNSTRANDED && c == '.')
			{
				PI p(ht->hid, hx->hid);
				bmap.insert(pair<PI, bam1_t>(p, b1t));

				bam1_t b1;
				bam1_t b2;
				build_bam1_t(b1, decode_vlist(ht->vlist), *ht, ht->pos, ht->rpos, ht->xs);
				build_bam1_t(b2, decode_vlist(hx->vlist), *hx, hx->pos, hx->rpos, hx->xs);
				umap.insert(pair<int, bam1_t>(ht->hid, b1));
				umap.insert(pair<int, bam1_t>(hx->hid, b2));
			}
			else
			{
				bam_write1(fout, &(b1t));
				assert(b1t.data != NULL);
				delete b1t.data;
			}

			ht = ht->next;
			hx = hx->next;
		}
	}

	for(int k = 0; k < hits.size(); k++)
	{
		hit &h = hits[k];
		if(h.bridged == true) continue;

		bam1_t b1t;
		build_bam1_t(b1t, decode_vlist(h.vlist), h, h.pos, h.rpos, h.xs);

		if(library_type == UNSTRANDED && h.xs == '.')
		{
			umap.insert(pair<int, bam1_t>(h.hid, b1t));
		}
		else
		{
			bam_write1(fout, &(b1t));
			assert(b1t.data != NULL);
			delete b1t.data;
		}

		/*
		if(h.nh >= 2)
		{
			uset.push_back(b1t);
		}
		else
		{
			bam_write1(fout, &(b1t));
			assert(b1t.data != NULL);
			delete b1t.data;
		}
		*/
	}
	return 0;
}
