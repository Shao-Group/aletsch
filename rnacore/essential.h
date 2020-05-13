/*
Part of meta-scallop
(c) 2020 by Mingfu Shao, The Pennsylvania State University
See LICENSE for licensing.
*/

#ifndef __ESSENTIAL_H__
#define __ESSENTIAL_H__

#include "hit.h"
#include "splice_graph.h"
#include "constants.h"

using namespace std;

// determine whether a->b is continuous in graph
bool check_continuous_vertices(splice_graph &gr, int x, int y);

// determine whether vv is a valid path in graph
bool check_valid_path(splice_graph &gr, const vector<int> &vv);

// construct sub-splice-graph
int build_child_splice_graph(splice_graph &root, splice_graph &gr, map<int, int> &ss);

// calculate total length of intron chains
int32_t get_total_length_of_introns(const vector<int32_t> &chain);

// transform between paths and coordinates
int build_exon_coordinates_from_path(splice_graph &gr, const vector<int> &v, vector<int32_t> &vv);
int build_intron_coordinates_from_path(splice_graph &gr, const vector<int> &v, vector<int32_t> &vv);
int check_strand_from_intron_coordinates(splice_graph &gr, const vector<int32_t> &v);
bool build_path_from_exon_coordinates(splice_graph &gr, const vector<int32_t> &v, vector<int> &vv);
bool build_path_from_intron_coordinates(splice_graph &gr, const vector<int32_t> &v, vector<int> &vv);
bool build_path_from_mixed_coordinates(splice_graph &gr, const vector<int32_t> &v, vector<int> &vv);

// align hit to splice graph
bool align_hit_to_splice_graph(const hit &h, splice_graph &gr, vector<int> &vv);

// match paired-end reads
int build_paired_reads(const vector<hit> &hits, vector<PI> &fs);

bool merge_intron_chains(const vector<int32_t> &x, const vector<int32_t> &y, vector<int32_t> &xy);
bool consistent_intron_chains(const vector<int32_t> &x, const vector<int32_t> &y);

// write bam
int add_cigar_skip(bam1_t &b1t, int32_t p1, int32_t p2);
int add_cigar_match(bam1_t &b1t, int32_t p1, int32_t p2);
int build_bam1_t(bam1_t &b1t, const hit &h);
int build_bam1_t(bam1_t &b1t, const hit &h1, const hit &h2, const vector<int32_t> &chain);

#endif
