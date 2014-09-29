
#include "util.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>

#include "codebook.h"
#include "qv_compressor.h"
#include "cluster.h"

/**
 *
 */
void encode(char *input_name, char *output_name, struct qv_options_t *opts) {
	struct quality_file_t qv_info;
	struct distortion_t *dist = generate_distortion_matrix(41, DISTORTION_MSE);
	struct alphabet_t *alphabet = alloc_alphabet(41);
	uint32_t status;
	struct hrtimer_t cluster_time, stats, encoding, total;
	FILE *fin, *fout;
	uint64_t bytes_used;
    double distortion;
	char tmp[100];

	start_timer(&total);
    
	// Load input file all at once
	status = load_file(input_name, &qv_info, 0);
	if (status != LF_ERROR_NONE) {
		printf("load_file returned error: %d\n", status);
		exit(1);
	}

	// Set up clustering data structures
	qv_info.alphabet = alphabet;
	qv_info.dist = dist;
	qv_info.cluster_count = opts->clusters;
	qv_info.clusters = alloc_cluster_list(&qv_info);
	qv_info.opts = opts;

scanf("%d", tmp);

	// Do k-means clustering
	start_timer(&cluster_time);
	do_kmeans_clustering(&qv_info);
	stop_timer(&cluster_time);
	if (opts->verbose) {
		printf("Clustering took %.4f seconds\n", get_timer_interval(&cluster_time));
	}
    
scanf("%d", tmp);

	// Then find stats and generate codebooks for each cluster
	start_timer(&stats);
	calculate_statistics(&qv_info);
	generate_codebooks(&qv_info);
	stop_timer(&stats);
    
scanf("%d", tmp);

	if (opts->verbose) {
		printf("Stats and codebook generation took %.4f seconds\n", get_timer_interval(&stats));
		// @todo expected distortion is inaccurate due to lack of pmf
		//printf("Expected distortion: %f\n", opts->e_dist);
	}
    
	// Note that we want \r\n translation in the input
	// but we do not want it in the output
	fin = fopen(input_name, "rt");
	fout = fopen(output_name, "wb");
	if (!fin || !fout) {
		perror("Unable to open input and output files");
		exit(1);
	}
	
	// @todo qv_compression should use quality_file structure with data in memory, now
	start_timer(&encoding);
	write_codebooks(fout, &qv_info);
    bytes_used = start_qv_compression(&qv_info, fout, &distortion);
	stop_timer(&encoding);
	stop_timer(&total);

	fclose(fin);
	fclose(fout);
    
	// Verbose stats
	if (opts->verbose) {
		// @todo add cluster info here
		printf("Actual distortion: %f\n", distortion);
		printf("Lines: %lu\n", qv_info.lines);
		printf("Columns: %u\n", qv_info.columns);
		printf("Total bytes used: %lu\n", bytes_used);
		printf("Encoding took %.4f seconds.\n", get_timer_interval(&total));
		printf("Total time elapsed: %.4f seconds.\n", get_timer_interval(&total));
	}

	// Parse-able stats
	if (opts->stats) {
		printf("rate, %.4f, distortion, %.4f, time, %.4f, size, %lu \n", (bytes_used*8.)/((double)(qv_info.lines)*qv_info.columns), distortion, get_timer_interval(&total), bytes_used);
	}
}

/**
 *
 */
void decode(char *input_file, char *output_file, struct qv_options_t *opts) {
	FILE *fin, *fout;
	struct hrtimer_t timer;
	struct quality_file_t qv_info;
	struct alphabet_t *A = alloc_alphabet(41);
    
	qv_info.alphabet = A;

	start_timer(&timer);

	fin = fopen(input_file, "rb");
	fout = fopen(output_file, "wt");
	if (!fin || !fout) {
		perror("Unable to open input or output files");
		exit(1);
	}

	read_codebooks(fin, &qv_info);
    start_qv_decompression(fout, fin, &qv_info);

	fclose(fout);
	fclose(fin);
	stop_timer(&timer);

	if (opts->verbose) {
		printf("Decoded %lu lines in %f seconds.\n", qv_info.lines, get_timer_interval(&timer));
	}
}

/**
 * Displays a usage name
 * @param name Program name string
 */
void usage(char *name) {
	printf("Usage: %s (options) [input file] [output file]\n", name);
	printf("Options are:\n");
	printf("\t-q\t\t: Store quality values in compressed file (default)\n");
	printf("\t-x\t\t: Extract quality values from compressed file\n");
	printf("\t-f [ratio]\t: Compress using [ratio] bits per bit of input entropy per symbol\n");
	printf("\t-r [rate]\t: Compress using fixed [rate] bits per symbol\n");
	printf("\t-c [#]\t: Compress using [#] clusters (default: 3)\n");
	printf("\t-h\t\t: Print this help\n");
	printf("\t-s\t\t: Print summary stats\n");
	printf("\t-t [lines]\t: Number of lines to use as training set (0 for all, 1000000 default)\n");
	printf("\t-v\t\t: Enable verbose output\n");
}

/**
 *
 */
int main(int argc, char **argv) {
    char *input_name = 0;
	char *output_name = 0;
	struct qv_options_t opts;
	uint32_t i;

	uint8_t extract = 0;
	uint8_t file_idx = 0;

	opts.training_size = 1000000;
	opts.verbose = 0;
	opts.stats = 0;
	opts.ratio = 0.5;
	opts.clusters = 3;

	// No dependency, cross-platform command line parsing means no getopt
	// So we need to settle for less than optimal flexibility (no combining short opts, maybe that will be added later)
	i = 1;
	while (i < argc) {
		// Handle file names and reject any other untagged arguments
		if (argv[i][0] != '-') {
			switch (file_idx) {
				case 0:
					input_name = argv[i];
					file_idx = 1;
					break;
				case 1:
					output_name = argv[i];
					file_idx = 2;
					break;
				default:
					printf("Garbage argument \"%s\" detected.\n", argv[i]);
					usage(argv[0]);
					exit(1);
			}
			i += 1;
			continue;
		}

		// Flags for options
		switch(argv[i][1]) {
			case 'x':
				extract = 1;
				i += 1;
				break;
			case 'q':
				extract = 0;
				i += 1;
				break;
			case 'f':
				extract = 0;
				opts.ratio = atof(argv[i+1]);
				opts.mode = MODE_RATIO;
				i += 2;
				break;
			case 'r':
				extract = 0;
				opts.ratio = atof(argv[i+1]);
				opts.mode = MODE_FIXED;
				i += 2;
				printf("--Warning-- fixed rate encoding not yet implemented, falling back to ratio");
				break;
			case 'c':
				opts.clusters = atoi(argv[i+1]);
				i += 2;
				break;
			case 'v':
				opts.verbose = 1;
				i += 1;
				break;
			case 'h':
				usage(argv[0]);
				exit(0);
			case 's':
				opts.stats = 1;
				i += 1;
				break;
			case 't':
				opts.training_size = atoi(argv[i+1]);
				i += 2;
				break;
			default:
				printf("Unrecognized option -%c.\n", argv[i][1]);
				usage(argv[0]);
				exit(1);
		}
	}

	if (file_idx != 2) {
		printf("Missing required filenames.\n");
		usage(argv[0]);
		exit(1);
	}

	if (opts.verbose) {
		if (extract) {
			printf("%s will be decoded to %s.\n", input_name, output_name);
		}
		else {
			printf("%s will be encoded as %s.\n", input_name, output_name);
			if (opts.mode == MODE_RATIO)
				printf("Ratio mode selected, targeting %f compression ratio\n", opts.ratio);
			else if (opts.mode == MODE_FIXED)
				printf("Fixed-rate mode selected, targeting %f bits per symbol\n", opts.ratio);
			else if (opts.mode == MODE_FIXED_MSE)
				printf("Fixed-MSE mode selected, targeting %f average MSE per context\n", opts.ratio);

			printf("Compression will use %d clusters\n", opts.clusters);
			// @todo other modes?
		}
	}

	if (extract) {
		decode(input_name, output_name, &opts);
	}
	else {
		encode(input_name, output_name, &opts);
	}

#ifndef LINUX
	system("pause");
#endif

	return 0;
}

