/*
Part of aletsch
(c) 2020 by Mingfu Shao, The Pennsylvania State University
See LICENSE for licensing.
*/

#ifndef __ESSENTIAL_H__
#define __ESSENTIAL_H__

#include "hit.h"
#include "hit_core.h"
#include "splice_graph.h"
#include "constants.h"
#include "pereads_cluster.h"
#include "phase_set.h"
#include "interval_map.h"

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
int annotate_path(splice_graph &gr, const vector<int32_t> &v, vector<int32_t> &vv, vector<int> &nn, vector<int> &pp);
int annotate_segment(splice_graph &gr, int32_t p1, int32_t p2, vector<int32_t> &vv, vector<int> &nn, vector<int> &pp);
int annotate_junction(splice_graph &gr, int32_t p1, int32_t p2, vector<int32_t> &vv, vector<int> &nn, vector<int> &pp);
bool build_path_from_exon_coordinates(splice_graph &gr, const vector<int32_t> &v, vector<int> &vv);
bool build_path_from_intron_coordinates(splice_graph &gr, const vector<int32_t> &v, vector<int> &vv);
bool build_path_from_mixed_coordinates(splice_graph &gr, const vector<int32_t> &v, vector<int> &vv);

// align hit to splice graph
bool align_hit_to_splice_graph(const hit &h, const vector<int32_t> &chain, splice_graph &gr, vector<int> &vv);

// check consistency
bool merge_intron_chains(const vector<int32_t> &x, const vector<int32_t> &y, vector<int32_t> &xy);
bool consistent_intron_chains(const vector<int32_t> &x, const vector<int32_t> &y);

// phase-set
int build_phase_set_from_unpaired_reads(phase_set &ps, splice_graph &gr, const vector<hit> &hits, const vector<bool> &paired);

// write bam
int add_cigar_skip(bam1_t &b1t, int32_t p1, int32_t p2);
int add_cigar_match(bam1_t &b1t, int32_t p1, int32_t p2);
bool build_bam1_t(bam1_t &b1t, const hit_core &h, const vector<int32_t> &chain);
bool build_bam1_t(bam1_t &b1t, const hit_core &h1, const hit_core &h2, const vector<int32_t> &chain);
int write_bridged_pereads_cluster(BGZF *fout, const pereads_cluster &pc, const vector<int32_t> &whole);
int write_unbridged_pereads_cluster(BGZF *fout, const pereads_cluster &pc);
int write_unpaired_reads(BGZF *fout, const vector<hit> &hits, const vector<bool> &paired);

// build transcript(s)
int build_transcript(splice_graph &gr, transcript &trst, const vector<int> &v, char strand, double abd, const string &tid, int type);
bool build_single_exon_transcript(splice_graph &gr, transcript &trst);

#endif
