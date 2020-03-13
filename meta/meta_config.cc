/*
Part of meta-scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#include "meta_config.h"
#include <cstdlib>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <cstring>

using namespace std;

// for meta-assembly
int min_supporting_samples = 2;			// 2
int min_splicing_count = 5;
int min_phasing_count = 1;
int32_t max_group_boundary_distance = 10000;
double max_group_junction_distance = 100;
bool merge_intersection = false;
int max_threads = 10;
int max_combined = 100;
double merge_threshold = 0.5;

int meta_verbose = 0;
string meta_version = "0.1.0";

string input_bam_list = "";
string output_gtf_file = "";

int parse_meta_arguments(int argc, const char ** argv)
{
	for(int i = 1; i < argc; i++)
	{
		// user specified

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
			printf("%s\n", meta_version.c_str());
			exit(0);
		}
		else if(string(argv[i]) == "--help")
		{
			print_meta_copyright();
			print_meta_help();
			printf("\n");
			print_meta_logo();
			exit(0);
		}
		else if(string(argv[i]) == "--merge_threshold")
		{
			merge_threshold = atof(argv[i + 1]);
			i++;
		}

		else if(string(argv[i]) == "--merge_intersection")
		{
			merge_intersection = true;
		}
		else if(string(argv[i]) == "--verbose")
		{
			meta_verbose = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--max_threads")
		{
			max_threads = atoi(argv[i + 1]);
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
		else if(string(argv[i]) == "--min_phasing_count")
		{
			min_phasing_count = atoi(argv[i + 1]);
			i++;
		}
	}
	return 0;
}

int print_meta_logo()
{
	return 0;
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

int print_meta_help()
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

int print_meta_copyright()
{
	printf("meta-scallop %s (c) 2020 Mingfu Shao, The Pennsylvania State University\n", meta_version.c_str());
	return 0;
}
