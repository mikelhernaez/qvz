
#include "util.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>

#include "codebook.h"
#include "qv_compressor.h"

/**
 *
 */
void encode(char *input_name, char *output_name, char *codebook_name, struct qv_options_t *opts) {
	// New variables, definitely used
	struct distortion_t *dist = generate_distortion_matrix(41, DISTORTION_MSE);
	struct alphabet_t *alphabet = alloc_alphabet(41);
	struct quality_file_t training_file;
	struct cond_pmf_list_t *training_stats;
	struct cond_quantizer_list_t *qlist;
	uint32_t status;
	struct hrtimer_t stats, encoding, total;
	FILE *fp;
	uint64_t bytes_used;
    double distortion;

	start_timer(&total);
	start_timer(&stats);
    
	// Load file statistics and find statistics for the training data
	status = load_file(input_name, &training_file, 1000000);
	if (status != LF_ERROR_NONE) {
		printf("load_file returned error: %d\n", status);
		exit(1);
	}
    
	training_stats = alloc_conditional_pmf_list(alphabet, training_file.columns);
	qlist = generate_codebooks(&training_file, training_stats, dist, opts);
    
	stop_timer(&stats);
	start_timer(&encoding);
    
	if (opts->verbose) {
		printf("Stats and codebook generation took %.4f seconds\n", get_timer_interval(&stats));
		printf("Expected distortion: %f\n", opts->e_dist);
	}
    
	fp = fopen(input_name, "rt");
	if (!fp) {
		perror("Unable to open input file");
		exit(1);
	}
	
    bytes_used = start_qv_compression(fp, output_name, qlist, &distortion);
	fclose(fp);

	// Write codebook second so that we know how many lines were encoded
	write_codebook(codebook_name, qlist);
	
	stop_timer(&encoding);
	stop_timer(&total);
    
	// Verbose stats
	if (opts->verbose) {
		printf("Actual distortion: %f\n", distortion);
		printf("Lines: %d\n", qlist->lines);
		printf("Total bytes used: %llu\n", bytes_used);
		printf("Encoding took %.4f seconds.\n", get_timer_interval(&total));
		printf("Total time elapsed: %.4f seconds.\n", get_timer_interval(&total));
	}

	// Parse-able stats
	if (opts->stats) {
		printf("rate, %.4f, distortion, %.4f, time, %.4f\n", (bytes_used*8.)/(qlist->lines*qlist->columns), distortion, get_timer_interval(&total));
	}
}

/**
 *
 */
void decode(char *input_file, char *output_file, char* codebook_file, struct qv_options_t *opts) {
	FILE *fout;
	struct hrtimer_t timer;
	struct cond_quantizer_list_t *qlist;
	struct alphabet_t *A = alloc_alphabet(41);
    
	start_timer(&timer);
    
	qlist = read_codebook(codebook_file, A);
	qlist->options = opts;
    
	fout = fopen(output_file, "wt");
	if (!fout) {
		perror("Unable to open output file");
		exit(1);
	}
    
    start_qv_decompression(fout, input_file, qlist);

	fclose(fout);
	stop_timer(&timer);

	if (opts->verbose) {
		printf("Decoded %d lines in %f seconds.\n", qlist->lines, get_timer_interval(&timer));
	}
}

/**
 * Displays a usage name
 * @param name Program name string
 */
void usage(char *name) {
	printf("Usage: %s (options) [codebook file] [input file] [output file]\n", name);
	printf("Options are:\n");
	printf("\t-x\t\t: Extract quality values from compressed file\n");
	printf("\t-c\t\t: Store quality values in compressed file (default)\n");
	printf("\t-f [ratio]\t: Compress using [ratio] bits per bit of input entropy per symbol\n");
	printf("\t-r [rate]\t: Compress using fixed [rate] bits per symbol\n");
	printf("\t-v\t\t: Enable verbose output\n");
	printf("\t-s\t\t: Print summary stats\n");
	printf("\t-h\t\t: Print this help\n");
}

/**
 *
 */
int main(int argc, char **argv) {
    char *input_name = 0;
	char *output_name = 0;
	char *codebook_name = 0;
	struct qv_options_t opts;
	uint32_t i;

	uint8_t extract = 0;
	uint8_t file_idx = 0;

	opts.verbose = 0;
	opts.stats = 0;
	opts.ratio = 0.5;

	// No dependency, cross-platform command line parsing means no getopt
	// So we need to settle for less than optimal flexibility (no combining short opts, maybe that will be added later)
	i = 1;
	while (i < argc) {
		// Handle file names and reject any other untagged arguments
		if (argv[i][0] != '-') {
			switch (file_idx) {
				case 0:
					codebook_name = argv[i];
					file_idx = 1;
					break;
				case 1:
					input_name = argv[i];
					file_idx = 2;
					break;
				case 2:
					output_name = argv[i];
					file_idx = 3;
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
			case 'c':
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
			default:
				printf("Unrecognized option -%c.\n", argv[i][1]);
				usage(argv[0]);
				exit(1);
		}
	}

	if (file_idx != 3) {
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
			// @todo other modes?
		}
	}

	if (extract) {
		decode(input_name, output_name, codebook_name, &opts);
	}
	else {
		encode(input_name, output_name, codebook_name, &opts);
	}

#ifndef LINUX
	system("pause");
#endif

	return 0;
}

