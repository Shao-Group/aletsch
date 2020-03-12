/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#ifndef __CONFIG_H__
#define __CONFIG_H__

#include "util.h"
#include <stdint.h>
#include <map>
#include <sstream>

using namespace std;

class config
{
public:
	config();

public:
	// for bam file and reads
	int min_flank_length;
	int max_num_cigar;
	int max_edit_distance;
	int32_t min_bundle_gap;
	int min_num_hits_in_bundle;
	uint32_t min_mapping_quality;
	bool uniquely_mapped_only;
	bool use_second_alignment;
	int min_splice_boundary_hits;

	// for preview
	bool preview_only;
	int max_preview_reads;
	int max_preview_spliced_reads;
	int min_preview_spliced_reads;
	double preview_infer_ratio;
	bool merge_intersection;

	// for identifying subgraphs
	int32_t min_subregion_gap;
	double min_subregion_overlap;
	int32_t min_subregion_length;
	int min_subregion_ladders;
	int32_t max_cluster_boundary_distance;
	int32_t max_cluster_intron_distance;
	double min_cluster_single_exon_ratio;

	// for subsetsum and router
	int max_dp_table_size;
	int min_router_count;

	// for splice graph
	double max_intron_contamination_coverage;
	double min_surviving_edge_weight;
	double max_decompose_error_ratio[7];
	double min_transcript_numreads;
	double min_transcript_coverage;
	double min_single_exon_coverage;
	double min_transcript_coverage_ratio; 
	int min_transcript_length_base;
	int min_transcript_length_increase;
	int min_exon_length;
	int max_num_exons;

	// for simulation
	int simulation_num_vertices;
	int simulation_num_edges;
	int simulation_max_edge_weight;

	// input and output
	string algo;
	string input_file;
	string graph_file;
	string output_file;

	// for controling
	bool output_tex_files;
	string fixed_gene_name;
	int max_num_bundles;
	int library_type;
	int min_gtf_transcripts_num;
	int batch_bundle_size;
	int verbose;
	string version;

public:
	int print_command_line(int argc, const char ** argv);
	int parse_arguments(int argc, const char ** argv);
	int print_copyright();
	int print_logo();
	int print_help();
};

#endif
