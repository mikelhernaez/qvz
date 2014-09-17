
#include "util.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>

#include "codebook.h"
#include "qv_compressor.h"

// Default file names if none are provided on the command line, for testing in Visual Studio
// @todo remove from final version
char *default_input = "quality_lines.txt";
char *default_output = "c_bitpacked_quality_lossy.txt";
char *default_codebook = "codebook.txt";

/**
 *
 */
void encode(char *input_name, char *output_name, char *codebook_name) {
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
	double distortion = 0.0;
    
	start_timer(&total);
	start_timer(&stats);
    
	// Load file statistics and find statistics for the training data
	status = load_file(input_name, &training_file, 1000000);
	if (status != LF_ERROR_NONE) {
		printf("load_file returned error: %d\n", status);
		exit(1);
	}
    
	training_stats = alloc_conditional_pmf_list(alphabet, training_file.columns);
	qlist = generate_codebooks(&training_file, training_stats, dist, 0.5, &distortion);
    
	stop_timer(&stats);
	start_timer(&encoding);
    
	printf("Stats and codebook generation took %.4f seconds\n", get_timer_interval(&stats));
	printf("Expected distortion: %f\n", distortion);
    
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
    
	// Summary stats from measurement
	printf("Actual distortion: %f\n", distortion);
	printf("Lines: %d\n", qlist->lines);
	printf("Total bytes used: %llu\n", bytes_used);
	printf("Encoding took %.4f seconds.\n", get_timer_interval(&total));
	printf("Total time elapsed: %.4f seconds.\n", get_timer_interval(&total));
}

/**
 *
 */
void decode(char *input_file, char *output_file, char* codebook_file) {
	FILE *fout;
	struct hrtimer_t timer;
	struct cond_quantizer_list_t *qlist;
	struct alphabet_t *A = alloc_alphabet(41);
    
	start_timer(&timer);
    
	qlist = read_codebook(codebook_file, A);
    
	fout = fopen(output_file, "wt");
	if (!fout) {
		perror("Unable to open output file");
		exit(1);
	}
    
    start_qv_decompression(fout, input_file, qlist);

	fclose(fout);
	stop_timer(&timer);
	printf("Decoded %d lines in %f seconds.\n", qlist->lines, get_timer_interval(&timer));
    
}

/**
 * Displays a usage name
 * @param name Program name string
 */
void usage(char *name) {
	printf("Usage: %s [mode] [codebook file] [input file] [output file]\n", name);
	printf("Where [mode] is -c for compression and -d for decompression\n");
}

/**
 *
 */
int main(int argc, char **argv) {
    char *input_name, *output_name, *codebook_name;
    
	// Check for arguments and get file names
	if (argc != 5) {
		usage(argv[0]);
		
		// Attempted to run from command line but not all arguments are present
		if (argc > 1) {
			printf("Incorrect number of arguments supplied!\n");
			exit(1);
		}
		else {
			printf("Using default input file: %s\nDefault output file: %s.\nDefault codebook file:%s\n", default_input, default_output, default_codebook);
			input_name = default_input;
			output_name = default_output;
            codebook_name = default_codebook;
		}
	}
	else {
        codebook_name = argv[2];
		input_name = argv[3];
		output_name = argv[4];
	}
    
    if (strcmp(argv[1], "-c") == 0) {
        encode(input_name, output_name, codebook_name);
    }
    else if (strcmp(argv[1], "-d") == 0) {
        decode(input_name, output_name, codebook_name);
    }
    else {
		usage(argv[0]);
		printf("Missing required -c or -d!\n");
		exit(1);
	}

#ifndef LINUX
	system("pause");
#endif

	return 0;
}

