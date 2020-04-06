/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#ifndef __PARAMETERS_H__
#define __PARAMETERS_H__

#include "util.h"
#include <stdint.h>
#include <map>
#include <sstream>

using namespace std;

class parameters
{
public:
	parameters();

public:
	// meta
	int min_supporting_samples;
	int min_splicing_count;
	int min_phasing_count;
	int32_t max_group_boundary_distance;
	bool merge_intersection;
	int max_threads;
	int max_combined;
	double merge_threshold;

	// for bam file and reads
	int min_flank_length;
	int max_num_cigar;
	int32_t min_bundle_gap;
	int min_num_hits_in_bundle;
	uint32_t min_mapping_quality;
	bool uniquely_mapped_only;
	bool use_second_alignment;

	// for preview
	bool preview_only;
	int max_preview_reads;
	int max_preview_spliced_reads;
	int min_preview_spliced_reads;
	double preview_infer_ratio;

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

	// for splice graph
	double max_intron_contamination_coverage;
	double min_surviving_edge_weight;
	double max_decompose_error_ratio[7];
	double min_transcript_coverage;
	double min_single_exon_coverage;
	int min_transcript_length_base;
	int min_transcript_length_increase;
	int min_exon_length;
	int max_num_exons;

	// input and output
	string algo;
	string input_file;
	string graph_file;
	string output_file;

	// for controling
	string input_bam_list;
	string output_gtf_file;
	int verbose;
	string version;
	int batch_bundle_size;

public:
	int print_command_line(int argc, const char ** argv);
	int parse_arguments(int argc, const char ** argv);
	int print_copyright();
	int print_logo();
	int print_help();
};

#endif
