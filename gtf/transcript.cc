/*
Part of aletsch
(c) 2020 by  Mingfu Shao, The Pennsylvania State University.
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#include <cstdio>
#include <cassert>
#include <sstream>
#include <algorithm>
#include <map>

#include "transcript.h"
#include "util.h"

transcript::transcript()
{
}

transcript::transcript(const item &e)
{
	assign(e);
	exons.clear();
}

transcript::~transcript()
{
}

int transcript::assign(const item &e)
{
	//assert(e.feature == "transcript");
	seqname = e.seqname;
	source = e.source;
	feature = e.feature;
	gene_id = e.gene_id;
	transcript_id = e.transcript_id;
	transcript_type = e.transcript_type;
	gene_type = e.gene_type;
	start = e.start;
	end = e.end;
	strand = e.strand;
	frame = e.frame;
	score = e.score;
	coverage = e.coverage;
	RPKM = e.RPKM;
	FPKM = e.FPKM;
	TPM = e.TPM;
	return 0;
}

bool transcript::operator< (const transcript &t) const
{
	int b = seqname.compare(t.seqname);
	if(b < 0) return true;
	if(b > 0) return false;
	if(exons.size() == 0) return true;
	if(t.exons.size() == 0) return false;
	if(exons[0].first < t.exons[0].first) return true;
	else return false;
}

int transcript::clear()
{
	exons.clear();
	seqname = "";
	source = "";
	feature = "";
	gene_id = "";
	transcript_id = "";
	transcript_type = "";
	gene_type = "";
	start = 0;
	end = 0;
	strand = '.';
	frame = -1;
	coverage = 0;
	score = 0;
	RPKM = 0;
	TPM = 0;
    cov2 = 0;
    conf = 0;
    abd = 0;
    count1 = 0;
    count2 = 0;

	return 0;
}

int transcript::add_exon(int s, int t)
{
	exons.push_back(PI32(s, t));
	return 0;
}

int transcript::add_exon(const item &e)
{
	assert(e.transcript_id == transcript_id);
	add_exon(e.start, e.end);
	return 0;
}

int transcript::sort()
{
	std::sort(exons.begin(), exons.end());
	return 0;
}

int transcript::shrink()
{
	if(exons.size() == 0) return 0;
	vector<PI32> v;
	PI32 p = exons[0];
	for(int i = 1; i < exons.size(); i++)
	{
		PI32 q = exons[i];
		if(p.second == q.first)
		{
			p.second = q.second;
		}
		else
		{
			//assert(p.second < q.first);
			v.push_back(p);
			p = q;
		}
	}
	v.push_back(p);
	exons = v;
	return 0;
}

int transcript::assign_RPKM(double factor)
{
	RPKM = coverage * factor;
	return 0;
}

int transcript::length() const
{
	int s = 0;
	for(int i = 0; i < exons.size(); i++)
	{
		assert(exons[i].second > exons[i].first);
		s += exons[i].second - exons[i].first;
	}
	return s;
}

PI32 transcript::get_bounds() const
{
	if(exons.size() == 0) return PI32(-1, -1);
	int32_t p = exons[0].first;
	int32_t q = exons[exons.size() - 1].second;
	return PI32(p, q);
}

PI32 transcript::get_first_intron() const
{
	if(exons.size() <= 1) return PI32(-1, -1);
	int32_t p = exons[0].second;
	int32_t q = exons[1].first;
	return PI32(p, q);
}

vector<PI32> transcript::get_intron_chain() const
{
	vector<PI32> v;
	if(exons.size() <= 1) return v;

	int32_t p = exons[0].second;
	for(int k = 1; k < exons.size(); k++)
	{
		int32_t q = exons[k].first;
		v.push_back(PI32(p, q));
		p = exons[k].second;
	}
	return v;
}

size_t transcript::get_intron_chain_hashing() const
{
	if(exons.size() == 0) return 0;

	if(exons.size() == 1)
	{
		size_t p = (exons[0].first + exons[0].second) / 10000;
		return p + 1;
	}

	vector<int32_t> vv;
	vector<PI32> v = get_intron_chain();
	for(int i = 0; i < v.size(); i++) 
	{
		vv.push_back(v[i].first);
		vv.push_back(v[i].second);
	}
	return vector_hash(vv) + 1;
}

bool transcript::intron_chain_match(const transcript &t) const
{
	if(exons.size() != t.exons.size()) return false;
	if(exons.size() <= 1) return false;
	int n = exons.size() - 1;
	if(exons[0].second != t.exons[0].second) return false;
	if(exons[n].first != t.exons[n].first) return false;
	for(int k = 1; k < n - 1; k++)
	{
		if(exons[k].first != t.exons[k].first) return false;
		if(exons[k].second != t.exons[k].second) return false;
	}
	return true;
}

int transcript::intron_chain_compare(const transcript &t) const
{
	if(exons.size() < t.exons.size()) return +1;
	if(exons.size() > t.exons.size()) return -1;
	if(exons.size() <= 1) return 0;

	int n = exons.size() - 1;
	if(exons[0].second < t.exons[0].second) return +1;
	if(exons[0].second > t.exons[0].second) return -1;
	for(int k = 1; k < n - 1; k++)
	{
		if(exons[k].first < t.exons[k].first) return +1;
		if(exons[k].first > t.exons[k].first) return -1;
		if(exons[k].second < t.exons[k].second) return +1;
		if(exons[k].second > t.exons[k].second) return -1;
	}
	if(exons[n].first < t.exons[n].first) return +1;
	if(exons[n].first > t.exons[n].first) return -1;
	return 0;
}

bool transcript::equal1(const transcript &t, double single_exon_overlap) const
{
	if(exons.size() != t.exons.size()) return false;

	if(seqname != t.seqname) return false;
	if(strand == '+' && t.strand == '-') return false;
	if(strand == '-' && t.strand == '+') return false;

	if(exons.size() == 1)
	{
		int32_t p1 = exons[0].first < t.exons[0].first ? exons[0].first : t.exons[0].first;
		int32_t p2 = exons[0].first < t.exons[0].first ? t.exons[0].first : exons[0].first;
		int32_t q1 = exons[0].second > t.exons[0].second ? exons[0].second : t.exons[0].second;
		int32_t q2 = exons[0].second > t.exons[0].second ? t.exons[0].second : exons[0].second;

		int32_t overlap = q2 - p2;
		if(overlap >= single_exon_overlap * length()) return true;
		if(overlap >= single_exon_overlap * t.length()) return true;
		return false;

		/*
		double overlap = (q2 - p2) * 1.0 / (q1 - p1);
		if(overlap < 0.8) return false;
		else return true;
		*/
	}

	return intron_chain_match(t);
}

int transcript::compare1(const transcript &t, double single_exon_overlap) const
{
	if(exons.size() < t.exons.size()) return +1;
	if(exons.size() > t.exons.size()) return -1;

	if(seqname < t.seqname) return +1;
	if(seqname > t.seqname) return -1;
	if(strand < t.strand) return +1;
	if(strand > t.strand) return -1;

	if(exons.size() == 1)
	{
		int32_t p1 = exons[0].first < t.exons[0].first ? exons[0].first : t.exons[0].first;
		int32_t p2 = exons[0].first < t.exons[0].first ? t.exons[0].first : exons[0].first;
		int32_t q1 = exons[0].second > t.exons[0].second ? exons[0].second : t.exons[0].second;
		int32_t q2 = exons[0].second > t.exons[0].second ? t.exons[0].second : exons[0].second;

		int32_t overlap = q2 - p2;
		if(overlap >= single_exon_overlap * length()) return 0;
		if(overlap >= single_exon_overlap * t.length()) return 0;

		//double overlap = (q2 - p2) * 1.0 / (q1 - p1);
		//if(overlap >= 0.8) return 0;

		if(exons[0].first < t.exons[0].first) return +1;
		if(exons[0].first > t.exons[0].first) return -1;
		if(exons[0].second < t.exons[0].second) return +1;
		if(exons[0].second > t.exons[0].second) return -1;
	}

	return intron_chain_compare(t);
}

int transcript::extend_bounds(const transcript &t)
{
	if(exons.size() == 0) return 0;
	if(t.exons.front().first < exons.front().first) exons.front().first = t.exons.front().first;
	if(t.exons.back().second > exons.back().second) exons.back().second = t.exons.back().second;
	return 0;
}

string transcript::label() const
{
	char buf[10240];
	PI32 p = get_bounds();
	sprintf(buf, "%s:%d-%d", seqname.c_str(), p.first, p.second);
	return string(buf);
}

int transcript::write(ostream &fout, double cov2, int count) const
{
	fout.precision(4);
	fout<<fixed;

	if(exons.size() == 0) return 0;
	
	PI32 p = get_bounds();

	fout<<seqname.c_str()<<"\t";				// chromosome name
	fout<<source.c_str()<<"\t";					// source
	fout<<"transcript\t";						// feature
	fout<<p.first + 1<<"\t";					// left position
	fout<<p.second<<"\t";						// right position
	fout<<1000<<"\t";							// score, now as expression
	fout<<strand<<"\t";							// strand
	fout<<".\t";								// frame
	fout<<"gene_id \""<<gene_id.c_str()<<"\"; ";
	fout<<"transcript_id \""<<transcript_id.c_str()<<"\"; ";
	if(gene_type != "") fout<<"gene_type \""<<gene_type.c_str()<<"\"; ";
	if(transcript_type != "") fout<<"transcript_type \""<<transcript_type.c_str()<<"\"; ";
	//fout<<"RPKM \""<<RPKM<<"\"; ";
	fout<<"cov \""<<coverage<<"\"; ";
	if(cov2 >= -0.5) fout<<"cov2 \""<<cov2<<"\"; ";
	if(count >= -0.5) fout<<"count \""<<count<<"\"; ";
	fout << endl;

	for(int k = 0; k < exons.size(); k++)
	{
		fout<<seqname.c_str()<<"\t";		// chromosome name
		fout<<source.c_str()<<"\t";			// source
		fout<<"exon\t";						// feature
		fout<<exons[k].first + 1<<"\t";		// left position
		fout<<exons[k].second<<"\t";		// right position
		fout<<1000<<"\t";					// score, now as expression
		fout<<strand<<"\t";					// strand
		fout<<".\t";						// frame
		fout<<"gene_id \""<<gene_id.c_str()<<"\"; ";
		fout<<"transcript_id \""<<transcript_id.c_str()<<"\"; ";
		fout<<"exon \""<<k + 1<<"\"; "<<endl;
	}
	return 0;
}

int transcript::write_features(ostream &stat_file) const
{
	/*
    ofstream stat_file;
    string filename;
    if(sample_id < 0) filename = "meta.trstFeature.csv";
    else filename = "gtf/"+ to_string(sample_id) + ".trstFeature.csv";
    stat_file.open(filename, fstream::app);
    stat_file.setf(ios::fixed, ios::floatfield);
    stat_file.precision(2);
	*/

    stat_file << transcript_id << '\t'       // Transcript ID
        << meta_tid << '\t'         //Transcript ID in meta.gtf
        << seqname << '\t'
        << coverage << '\t'        // Meta coverage
        << cov2 << '\t'            // Individual coverage 2 
        << abd << '\t'             // Abundance
        << conf << '\t'            // Confidence
        << count1 << '\t'          // Actual count of meta trst       
        << count2 << '\t'          // Inferred count by individual trst
        << exons.size() << '\t'   // Number of exons
        << features.gr_vertices << '\t'
        << features.gr_edges << '\t'
        << features.gr_reads << '\t'
        << features.gr_subgraph << '\t'
        << features.num_vertices << '\t'
        << features.num_edges << '\t'
        << features.junc_ratio << '\t'
        << features.max_mid_exon_len << '\t'
        << features.start_loss1 << '\t'
        << features.start_loss2 << '\t'
        << features.start_loss3 << '\t'
        << features.end_loss1 << '\t'
        << features.end_loss2 << '\t'
        << features.end_loss3 << '\t'
        << features.start_merged_loss << '\t'
        << features.end_merged_loss << '\t'
        << features.introns << '\t'
        << features.intron_ratio << '\t'
        << features.start_introns << '\t'
        << features.start_intron_ratio << '\t'
        << features.end_introns << '\t'
        << features.end_intron_ratio << '\t'
        << features.uni_junc << '\t'
        << features.seq_min_wt << '\t'
        << features.seq_min_cnt << '\t'
        << features.seq_min_abd << '\t'
        << features.seq_min_ratio << '\t'
        << features.seq_max_wt << '\t'
        << features.seq_max_cnt << '\t'
        << features.seq_max_abd << '\t'
        << features.seq_max_ratio << '\t'
        << features.start_cnt << '\t'
        << features.start_weight << '\t'
        << features.start_abd << '\t'
        << features.end_cnt << '\t'
        << features.end_weight << '\t'
        << features.end_abd << '\t'
        << features.unbridge_start_coming_count << '\t'
        << features.unbridge_start_coming_ratio << '\t'
        << features.unbridge_end_leaving_count << '\t'
        << features.unbridge_end_leaving_ratio << endl;
        
    //stat_file.close();
    return 0;
}

int transcript::write_features(int sample_id) const
{
    ofstream stat_file;
    string filename;
    if(sample_id < 0) filename = "meta.trstFeature.csv";
    else filename = "gtf/"+ to_string(sample_id) + ".trstFeature.csv";

    stat_file.open(filename, fstream::app);
    stat_file.setf(ios::fixed, ios::floatfield);
    stat_file.precision(2);
    stat_file << transcript_id << '\t'       // Transcript ID
        << meta_tid << '\t'         //Transcript ID in meta.gtf
        << seqname << '\t'
        << coverage << '\t'        // Meta coverage
        << cov2 << '\t'            // Individual coverage 2 
        << abd << '\t'             // Abundance
        << conf << '\t'            // Confidence
        << count1 << '\t'          // Actual count of meta trst       
        << count2 << '\t'          // Inferred count by individual trst
        << exons.size() << '\t'   // Number of exons
        << features.gr_vertices << '\t'
        << features.gr_edges << '\t'
        << features.gr_reads << '\t'
        << features.gr_subgraph << '\t'
        << features.num_vertices << '\t'
        << features.num_edges << '\t'
        << features.junc_ratio << '\t'
        << features.max_mid_exon_len << '\t'
        << features.start_loss1 << '\t'
        << features.start_loss2 << '\t'
        << features.start_loss3 << '\t'
        << features.end_loss1 << '\t'
        << features.end_loss2 << '\t'
        << features.end_loss3 << '\t'
        << features.start_merged_loss << '\t'
        << features.end_merged_loss << '\t'
        << features.introns << '\t'
        << features.intron_ratio << '\t'
        << features.start_introns << '\t'
        << features.start_intron_ratio << '\t'
        << features.end_introns << '\t'
        << features.end_intron_ratio << '\t'
        << features.uni_junc << '\t'
        << features.seq_min_wt << '\t'
        << features.seq_min_cnt << '\t'
        << features.seq_min_abd << '\t'
        << features.seq_min_ratio << '\t'
        << features.seq_max_wt << '\t'
        << features.seq_max_cnt << '\t'
        << features.seq_max_abd << '\t'
        << features.seq_max_ratio << '\t'
        << features.start_cnt << '\t'
        << features.start_weight << '\t'
        << features.start_abd << '\t'
        << features.end_cnt << '\t'
        << features.end_weight << '\t'
        << features.end_abd << '\t'
        << features.unbridge_start_coming_count << '\t'
        << features.unbridge_start_coming_ratio << '\t'
        << features.unbridge_end_leaving_count << '\t'
        << features.unbridge_end_leaving_ratio << endl;
        
    stat_file.close();
    return 0;
}

void transcript::write_seq_features(ofstream & stat_file, const vector<int>& v) const
{
    stat_file << '[';
    for(size_t i = 0; i < v.size(); i++)
    {
        stat_file << v[i];
        if(i < v.size()-1) stat_file << ',';
    }
    stat_file << "]\t";
}

void transcript::write_seq_features(ofstream & stat_file, const vector<double>& v) const
{
    stat_file << '[';
    for(size_t i = 0; i < v.size(); i++)
    {
        stat_file << v[i];
        if(i < v.size()-1) stat_file << ',';
    }
    stat_file << "]\t";
}
