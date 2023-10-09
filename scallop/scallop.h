/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#ifndef __SCALLOP3_H__
#define __SCALLOP3_H__

#include "splice_graph.h"
#include "hyper_set.h"
#include "equation.h"
#include "router.h"
#include "path.h"
#include "parameters.h"

typedef map<edge_descriptor, double> MED;
typedef map<edge_descriptor, int> MEI;
typedef map< edge_descriptor, vector<int> > MEV;
typedef pair< edge_descriptor, vector<int> > PEV;
typedef pair< edge_descriptor, double > PED;
typedef pair< edge_descriptor, int > PEI;
typedef pair< vector<int>, vector<int> > PVV;
typedef pair<PEE, int> PPEEI;
typedef map<PEE, int> MPEEI;
typedef pair<int, int> PI;
typedef map<int, int> MI;

// for noisy splice graph
class scallop
{
public:
	scallop(splice_graph &gr, hyper_set &hs, const parameters &c, bool random_ordering = false);
	virtual ~scallop();

public:
	int assemble();

public:
	const parameters &cfg;				// parameters
	bool random_ordering;				// if shuffling vertices
	splice_graph &gr;					// splice graph
	hyper_set &hs;						// hyper edges
	MEI e2i;							// edge map, from edge to index
	VE i2e;								// edge map, from index to edge
	MEV mev;							// super edges
	MED med;							// map from super edges to ave-reads
	MEI mei;							// map from super edges to length
	vector<int> v2v;					// vertex map
	int round;							// iteration
	set<int> nonzeroset;				// vertices with degree >= 1
	vector<path> paths;					// predicted paths
	vector<transcript> trsts;			// predicted transcripts

private:
	// init
	int classify();
	int init_vertex_map();
	int init_super_edges();
	int init_inner_weights();
	int init_nonzeroset();
	int add_pseudo_hyper_edges();
	int refine_splice_graph();

	// resolve iteratively
	bool resolve_broken_vertex();
	bool target_smallest_edges(double max_ratio);
	bool remove_smallest_edges(double max_ratio);
	bool thread_smallest_edges(double max_ratio, int degree);
	bool resolve_smallest_edges(double max_ratio);
	bool resolve_trivial_vertex(int type, bool fast, double jump_ratio);
	bool resolve_single_trivial_vertex(int i, double jump_ratio);
	bool resolve_splittable_vertex(int type, int degree, double max_ratio);
	bool resolve_unsplittable_vertex(int type, int degree, double max_ratio);
	bool resolve_hyper_edge(int fsize);
	bool resolve_mixed_vertex(int type);
	bool resolve_mixed_smallest_edges();
	bool resolve_trivial_vertex_fast(double jump_ratio);
	bool resolve_single_trivial_vertex_fast(int i, double jump_ratio);

	// smooth vertex
	int balance_vertex(int x);
	int balance_vertex(int v, const vector<int> &ve1, const vector<int> &ve2);
	double compute_balance_ratio(int x);

	// decomposing subroutines
	bool remove_single_smallest_edge(int root, double max_ratio, double &ratio);
	bool remove_single_smallest_in_edge(int root, double max_ratio, double &ratio);
	bool remove_single_smallest_out_edge(int root, double max_ratio, double &ratio);
	bool thread_single_smallest_edge(int root, double max_ratio, double &ratio);
	bool thread_single_smallest_edge(int root, double max_ratio, double &ratio, int degree);
	int target_single_smallest_in_edge(int root, double max_ratio);
	int target_single_smallest_out_edge(int root, double max_ratio);
	int compute_smallest_edge(int x, double &ratio);
	int compute_smallest_in_edge(int x, double &ratio);
	int compute_smallest_out_edge(int x, double &ratio);
	int compute_second_smallest_in_edge(int x, double &ratio);
	int compute_second_smallest_out_edge(int x, double &ratio);
    int compute_smallest_edge_sample_abundance(int x);
    bool continuous_vertices(edge_descriptor e);
    bool closed_vertex(edge_descriptor e, int root);
    int update_log_confidence(int root);
	int decompose_trivial_vertex(int v);
	int decompose_vertex_extend(int v, MPID &pe2w);
	int decompose_vertex_replace(int v, MPID &pe2w);
	int classify_trivial_vertex(int v, bool fast);
	int exchange_sink(int old_sink, int new_sink);
	int split_vertex(int x, const vector<int> &xe, const vector<int> &ye);
	int split_edge(int exi, double w);
	int merge_adjacent_edges(int x, int y);
	int merge_adjacent_edges(int x, int y, double ww);
	int merge_adjacent_equal_edges(int x, int y);
	int remove_edge(int e);
	int split_merge_path(const VE &p, double w);
	int split_merge_path(const vector<int> &p, double w);
	int terminate_blocked_edges(int root);
	int terminate_smallest_edge(int root);
	int terminate_edge_sink(edge_descriptor e);
	int terminate_edge_source(edge_descriptor e);
	int borrow_edge_strand(int e1, int e2);
	bool consistent_strands(int e1, int e2);
	int break_divided_phases();
	bool left_adjacent(edge_descriptor e);
	bool right_adjacent(edge_descriptor e);

	int collect_existing_st_paths();
	int collect_path(int e);
	int compute_length(const path &p);
	int greedy_decompose();

	// stats, print, and draw
	int print_super_edges();
	int print_phasing_paths(hyper_set &hh);
	int print();
	int stats();
	int summarize_vertices();
	int draw_splice_graph(const string &file);
	vector<int> topological_sort();

	// collect paths and build transcripts
	int collect_phasing_paths();
	int collect_phasing_path(int e, int s, int t);
	int build_transcripts(splice_graph &gr);
    int update_trst_features(splice_graph &gr, transcript &trst, int i, vector<path> &paths);
    int check_junc_relation(const vector<pair<int,int>>& junc1, const vector<pair<int,int>>& junc2);
};

#endif
