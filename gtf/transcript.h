/*
Part of aletsch
(c) 2020 by  Mingfu Shao, The Pennsylvania State University.
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#ifndef __GTF_TRANSCRIPT_H__
#define __GTF_TRANSCRIPT_H__

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <set>
#include "item.h"

using namespace std;

typedef pair<int32_t, int32_t> PI32;

class transcript
{
public:
	transcript(const item &ie);
	transcript();
	~transcript();

public:
	bool operator< (const transcript &t) const;

public:
	string seqname;
	string source;
	string feature;
	string gene_id;
	string transcript_id;
	string gene_type;
	string transcript_type;
	int32_t start;
	int32_t end;
	double score;
	char strand;
	int frame;
	double coverage;
	double covratio;
	double RPKM;
	double FPKM;
	double TPM;

    //other features
    string meta_tid;
    double cov2;//coverge in individual sample
    double conf;//confidence in individual sample
    double abd;//overall abundance inferred from individuals
    int count1;// inferred count of samples supporting the trst
    int count2;// actual count of trst in output

    // Store the features for individual trst
    struct TrstFeatures {
    int gr_vertices;     // Number of vertices in the graph
    int gr_edges;        // Number of edges in the graph
    int gr_reads;
    int gr_subgraph;
    int num_vertices;       // Number of vertices in the path
    int num_edges;      //Number of edges in the path
    double junc_ratio;       //Ratio of junctions(length>1) in the path, except starting and ending
    int max_mid_exon_len;
    double start_loss1;
    double start_loss2;
    double start_loss3;
    double end_loss1;
    double end_loss2;
    double end_loss3;
    double start_merged_loss;
    double end_merged_loss;
    int introns;
    int start_introns;
    int end_introns;
    double intron_ratio; 
    double start_intron_ratio;
    double end_intron_ratio;
    int uni_junc;
    double seq_min_wt;
    int seq_min_cnt;
    double seq_min_abd;
    double seq_min_ratio;
    double seq_max_wt;
    int seq_max_cnt;
    double seq_max_abd;
    double seq_max_ratio;
    int unbridge_start_coming_count;
    double unbridge_start_coming_ratio;
    int unbridge_end_leaving_count;
    double unbridge_end_leaving_ratio;
    int start_cnt;
    double start_weight;
    double start_abd;
    int end_cnt;
    double end_weight;
    double end_abd;

    };

    TrstFeatures features;

	vector<PI32> exons;

public:
	int add_exon(int s, int t);
	int add_exon(const item &e);
	int assign_RPKM(double factor);
	int sort();
	int clear();
	int shrink();
	int assign(const item &e);
	int length() const;
	PI32 get_bounds() const;
	PI32 get_first_intron() const;
	vector<PI32> get_intron_chain() const;
	size_t get_intron_chain_hashing() const;
	bool intron_chain_match(const transcript &t) const;
	int intron_chain_compare(const transcript &t) const;
	bool equal1(const transcript &t, double single_exon_overlap) const;
	int compare1(const transcript &t, double single_exon_overlap) const;
	int extend_bounds(const transcript &t);
	string label() const;
	int write(ostream &fout, double cov2 = -1, int count = -1) const;
    int write_features(int sample_id) const;
    void write_seq_features(ofstream & stat_file, const vector<int>& v) const;
    void write_seq_features(ofstream & stat_file, const vector<double>& v) const;

};

#endif
