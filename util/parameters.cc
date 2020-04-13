/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#include "parameters.h"
#include <cstdlib>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <cstring>

using namespace std;

parameters::parameters()
{
	// for controling
	input_bam_list = "";
	output_gtf_file = "";
	verbose = 1;
	algo = "meta-scallop";
	version = "0.1.1";
	max_threads = 10;

	// for meta-assembly
	min_supporting_samples = 2;	
	min_splicing_count = 5;
	min_phasing_count = 1;
	max_combined = 100;
	merge_threshold = 0.5;
	standalone_coverage = 10.0;

	// for bridging paired-end reads
	bridge_dp_solution_size = 10;
	bridge_dp_stack_size = 5;

	// for loading bam file and reads
	min_flank_length = 3;
	max_num_cigar = 1000;
	min_bundle_gap = 50;
	min_num_hits_in_bundle = 20;
	min_mapping_quality = 1;
	use_second_alignment = false;
	uniquely_mapped_only = false;
	batch_bundle_size = 100;
	max_reads_partition_gap = 10;
	
	// for preview
	max_preview_reads = 2000000;
	max_preview_spliced_reads = 50000;
	min_preview_spliced_reads = 10000;
	preview_infer_ratio = 0.8;
	
	// for identifying subgraphs
	min_subregion_gap = 3;
	min_subregion_overlap = 1.5;
	min_subregion_length = 15;
	
	// for revising splice graph and phasing paths
	min_confident_phasing_path_weight = 2;
	max_group_boundary_distance = 10000;
	max_intron_contamination_coverage = 2.0;
	min_surviving_edge_weight = 1.5;

	// for decomposing splice graph
	max_decompose_error_ratio[0] = 0.33;
	max_decompose_error_ratio[1] = 0.05;
	max_decompose_error_ratio[2] = 0.00;
	max_decompose_error_ratio[3] = 0.25;
	max_decompose_error_ratio[4] = 0.30;
	max_decompose_error_ratio[5] = 0.00;
	max_decompose_error_ratio[6] = 1.10;
	min_guaranteed_edge_weight = 0.01;
	max_dp_table_size = 10000;
	
	// for filtering paths
	min_transcript_coverage = 1.01;
	min_single_exon_coverage = 20;
	min_transcript_length_base = 150;
	min_transcript_length_increase = 50;
	min_exon_length = 20;
	max_num_exons = 10000;

	// for clustering assembled transcripts
	max_cluster_boundary_distance = 10000;
	max_cluster_intron_distance = 5;
	min_cluster_single_exon_ratio = 0.8;
}

int parameters::parse_arguments(int argc, const char ** argv)
{
	for(int i = 1; i < argc; i++)
	{
		if(string(argv[i]) == "-i")
		{
			input_bam_list = string(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "-o")
		{
			output_gtf_file = string(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "-t")
		{
			max_threads = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--version")
		{
			printf("%s\n", version.c_str());
			exit(0);
		}
		else if(string(argv[i]) == "--help")
		{
			print_copyright();
			print_help();
			printf("\n");
			print_logo();
			exit(0);
		}
		else if(string(argv[i]) == "--verbose")
		{
			verbose = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--merge_threshold")
		{
			merge_threshold = atof(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--standalone_coverage")
		{
			standalone_coverage = atof(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--min_supporting_samples")
		{
			min_supporting_samples = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--min_splicing_count")
		{
			min_splicing_count = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--max_combined")
		{
			max_combined = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--bridge_dp_solution_size")
		{
			bridge_dp_solution_size = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--bridge_dp_stack_size")
		{
			bridge_dp_stack_size = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--min_phasing_count")
		{
			min_phasing_count = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--max_reads_partition_gap")
		{
			max_reads_partition_gap = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--min_flank_length")
		{
			min_flank_length = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--max_num_cigar")
		{
			max_num_cigar = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--min_bundle_gap")
		{
			min_bundle_gap = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--min_num_hits_in_bundle")
		{
			min_num_hits_in_bundle = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--min_mapping_quality")
		{
			min_mapping_quality = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--max_preview_spliced_reads")
		{
			max_preview_spliced_reads = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--min_preview_spliced_reads")
		{
			min_preview_spliced_reads = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--max_preview_reads")
		{
			max_preview_reads = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--preview_infer_ratio")
		{
			preview_infer_ratio = atof(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--min_subregion_gap")
		{
			min_subregion_gap = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--min_subregion_length")
		{
			min_subregion_length = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--min_subregion_overlap")
		{
			min_subregion_overlap = atof(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--min_surviving_edge_weight")
		{
			min_surviving_edge_weight = atof(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--max_intron_contamination_coverage")
		{
			max_intron_contamination_coverage = atof(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--min_transcript_coverage")
		{
			min_transcript_coverage = atof(argv[i + 1]);
			i++;
			if(fabs(min_transcript_coverage - 1.0) < 0.01) min_transcript_coverage = 1.01;
		}
		else if(string(argv[i]) == "--min_single_exon_coverage")
		{
			min_single_exon_coverage = atof(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--min_transcript_length_base")
		{
			min_transcript_length_base = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--min_transcript_length_increase")
		{
			min_transcript_length_increase = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--min_exon_length")
		{
			min_exon_length = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--max_num_exons")
		{
			max_num_exons = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--max_dp_table_size")
		{
			max_dp_table_size = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--max_decompose_error_ratio0")
		{
			max_decompose_error_ratio[0] = atof(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--max_decompose_error_ratio1")
		{
			max_decompose_error_ratio[1] = atof(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--max_decompose_error_ratio2")
		{
			max_decompose_error_ratio[2] = atof(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--max_decompose_error_ratio3")
		{
			max_decompose_error_ratio[3] = atof(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--max_decompose_error_ratio4")
		{
			max_decompose_error_ratio[4] = atof(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--max_decompose_error_ratio5")
		{
			max_decompose_error_ratio[5] = atof(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--max_decompose_error_ratio6")
		{
			max_decompose_error_ratio[6] = atof(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--use_second_alignment")
		{
			string s(argv[i + 1]);
			if(s == "true") use_second_alignment = true;
			else use_second_alignment = false;
			i++;
		}
		else if(string(argv[i]) == "--uniquely_mapped_only")
		{
			string s(argv[i + 1]);
			if(s == "true") uniquely_mapped_only = true;
			else uniquely_mapped_only = false;
			i++;
		}
		else if(string(argv[i]) == "--verbose")
		{
			verbose = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--batch_bundle_size")
		{
			batch_bundle_size = atoi(argv[i + 1]);
			i++;
		}
	}

	if(min_surviving_edge_weight < 0.1 + min_transcript_coverage) 
	{
		min_surviving_edge_weight = 0.1 + min_transcript_coverage;
	}

	return 0;
}

int parameters::print_command_line(int argc, const char ** argv)
{
	printf("command line: ");
	for(int i = 0; i < argc; i++)
	{
		printf("%s ", argv[i]);
	}
	printf("\n");
	return 0;
}

int parameters::print_logo()
{
	printf("      ___           ___           ___                                       ___           ___    \n");
	printf("     /  /\\         /  /\\         /  /\\                                     /  /\\         /  /\\   \n");
	printf("    /  /:/_       /  /:/        /  /::\\                                   /  /::\\       /  /::\\  \n");
	printf("   /  /:/ /\\     /  /:/        /  /:/\\:\\    ___     ___   ___     ___    /  /:/\\:\\     /  /:/\\:\\ \n");
	printf("  /  /:/ /::\\   /  /:/  ___   /  /:/~/::\\  /__/\\   /  /\\ /__/\\   /  /\\  /  /:/  \\:\\   /  /:/~/:/ \n");
	printf(" /__/:/ /:/\\:\\ /__/:/  /  /\\ /__/:/ /:/\\:\\ \\  \\:\\ /  /:/ \\  \\:\\ /  /:/ /__/:/ \\__\\:\\ /__/:/ /:/  \n");
	printf(" \\  \\:\\/:/~/:/ \\  \\:\\ /  /:/ \\  \\:\\/:/__\\/  \\  \\:\\  /:/   \\  \\:\\  /:/  \\  \\:\\ /  /:/ \\  \\:\\/:/   \n");
	printf("  \\  \\::/ /:/   \\  \\:\\  /:/   \\  \\::/        \\  \\:\\/:/     \\  \\:\\/:/    \\  \\:\\  /:/   \\  \\::/    \n");
	printf("   \\__\\/ /:/     \\  \\:\\/:/     \\  \\:\\         \\  \\::/       \\  \\::/      \\  \\:\\/:/     \\  \\:\\    \n");
	printf("     /__/:/       \\  \\::/       \\  \\:\\         \\__\\/         \\__\\/        \\  \\::/       \\  \\:\\   \n");
	printf("     \\__\\/         \\__\\/         \\__\\/                                     \\__\\/         \\__\\/   \n");
	printf("\n");

	return 0;
}

int parameters::print_help()
{
	printf("\n");
	printf("Usage: meta-scallop -i <input-bam-list> -o <output.gtf> [options]\n");
	printf("\n");
	printf("Options:\n");
	printf(" %-42s  %s\n", "--help",  "print usage of meta-scallop and exit");
	printf(" %-42s  %s\n", "--version",  "print current version of meta-scallop and exit");
	printf(" %-42s  %s\n", "-t/--max_threads <integer>",  "number of threads, default 10");
	printf(" %-42s  %s\n", "--max_combined <integer>",  "the maximized number of splice graphs that will be combined, default: 100");
	printf(" %-42s  %s\n", "--merge_threshold <float>",  "the minimized similarity for two graphs to be combined, default: 0.5");
	printf(" %-42s  %s\n", "--min_supporting_samples <integer>",  "the minimized number of samples needed to support a splicing site, default: 2");
	printf(" %-42s  %s\n", "--min_splicing_count <integer>",  "the minimized coverage needed to support a splicing site, default: 5");
	printf(" %-42s  %s\n", "--min_splice_bundary_hits <integer>",  "the minimum number of spliced reads required to support a junction, default: 1");
	printf(" %-42s  %s\n", "--min_transcript_coverage <float>",  "minimum coverage required for a multi-exon transcript, default: 1.01");
	printf(" %-42s  %s\n", "--min_single_exon_coverage <float>",  "minimum coverage required for a single-exon transcript, default: 20");
	printf(" %-42s  %s\n", "--min_transcript_length_increase <integer>",  "default: 50");
	printf(" %-42s  %s\n", "--min_transcript_length_base <integer>",  "default: 150, minimum length of a transcript would be");
	printf(" %-42s  %s\n", "",  "--min_transcript_length_base + --min_transcript_length_increase * num-of-exons");
	printf(" %-42s  %s\n", "--min_mapping_quality <integer>",  "ignore reads with mapping quality less than this value, default: 1");
	printf(" %-42s  %s\n", "--max_num_cigar <integer>",  "ignore reads with CIGAR size larger than this value, default: 1000");
	printf(" %-42s  %s\n", "--min_bundle_gap <integer>",  "minimum distances required to start a new bundle, default: 50");
	printf(" %-42s  %s\n", "--min_num_hits_in_bundle <integer>",  "minimum number of reads required in a bundle, default: 20");
	printf(" %-42s  %s\n", "--min_flank_length <integer>",  "minimum match length in each side for a spliced read, default: 3");
	return 0;
}

int parameters::print_copyright()
{
	printf("meta-scallop %s (c) 2020 Mingfu Shao, The Pennsylvania State University\n", version.c_str());
	return 0;
}
