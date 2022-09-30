/*
Part of aletsch
(c) 2020 by Mingfu Shao, The Pennsylvania State University
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
	// for controling
	string input_bam_list;
	string chrm_list_file;
	string chrm_list_string;
	string output_gtf_file;
	string output_gtf_dir;
	string output_bridged_bam_dir;
	string profile_dir;
	int verbose;
	string algo;
	string version;
	int max_threads;
	bool profile_only;
	bool boost_precision;

	// for meta-assembly
	int max_group_size;
	double min_grouping_similarity;
	double max_grouping_similarity;
	int max_num_junctions_to_combine;
	int assembly_repeats;
	//int32_t region_partition_length;

	// for bridging paired-end reads
	int bridge_end_relaxing;
	int bridge_dp_solution_size;
	int bridge_dp_stack_size;
	double min_bridging_score;

	// for decomposing splice graph
	double max_decompose_error_ratio[8];
	double min_guaranteed_edge_weight;
	int max_dp_table_size;

	// for filtering paths
	double min_transcript_coverage;
	double min_single_exon_transcript_coverage;
	double min_single_exon_individual_coverage;
	int min_transcript_length_base;
	int min_transcript_length_increase;
	int min_single_exon_transcript_length;
	int min_exon_length;
	int max_num_exons;

	// for clustering assembled transcripts
	int32_t max_cluster_boundary_distance;
	int32_t max_cluster_intron_distance;
	int32_t max_cluster_intron_shifting;
	double long_reads_cluster_boosting;
	double min_single_exon_clustering_overlap;

	// for loading bam file and reads
	int min_flank_length;
	int max_num_cigar;
	int32_t min_bundle_gap;
	int min_num_hits_in_bundle;
	int min_junction_support;
	uint32_t min_mapping_quality;
	bool uniquely_mapped_only;
	bool use_second_alignment;
	int batch_bundle_size;
	int32_t max_reads_partition_gap;

	// for preview
	int max_preview_reads;
	int max_preview_spliced_reads;
	int min_preview_spliced_reads;
	double preview_infer_ratio;

	// for identifying subgraphs
	int32_t min_subregion_gap;
	int32_t min_subregion_length;
	double min_subregion_overlap;
	double min_subregion_pvalue;

	// for revising splice graph and phasing paths
	int32_t max_group_boundary_distance;
	double max_intron_contamination_coverage;
	double min_surviving_edge_weight;
	int normal_junction_threshold;
	int extend_junction_threshold;
	double min_boundary_log_ratio;

public:
	int print_command_line(int argc, const char ** argv);
	int parse_arguments(int argc, const char ** argv, int data_type);
	int print_copyright();
	int print_logo();
	int print_help();

	int set_default(int data_type);
};

#endif
