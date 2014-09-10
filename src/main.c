
#define _CRT_SECURE_NO_WARNINGS

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>

#ifdef LINUX
	#define _stat stat
#endif

#include "codebook.h"
#include "util.h"

#define SYMBOLS 41

// Default file names if none are provided on the command line
const char *default_input = "quality_lines.txt";
const char *default_output = "c_bitpacked_quality_lossy.txt";

/**
 * Actually does the encoding (and also tests decoding to ensure correctness)
 */
int main(int argc, char **argv) {
	FILE *fp, *fbitout; //, *fref;
	uint32_t columns;

	// Old variables, might not all be needed as we improve the implementation
	uint32_t bits_used;
	uint32_t j = 0, s = 0, u = 0, k = 0;
	double distortion = 0.0;
	uint32_t error = 0;
	uint8_t state_enc = 0;
	uint8_t mask = 0;
	char *line;
	uint8_t *outline, *stateline, *precodeline;
	uint8_t *bits;

	// New variables, definitely used
	struct distortion_t *dist = generate_distortion_matrix(41, DISTORTION_MSE);
	struct alphabet_t *alphabet = alloc_alphabet(41);
	struct quality_file_t training_file;
	struct cond_pmf_list_t *training_stats;
	struct cond_quantizer_list_t *qlist;
	uint32_t status;
	struct hrtimer_t stats, encoding, total;
	const char *input_name, *output_name;
	struct quantizer_t *q;

	// Check for arguments and get file names
	if (argc != 3) {
		// Print usage guidelines
		printf("Usage: %s [input file] [output file]\n", argv[0]);
		
		// Attempted to run from command line but not all arguments are present
		if (argc > 1) {
			printf("Incorrect number of arguments supplied!");
			exit(1);
		}
		else {
			printf("Using default input file: %s\nDefault output file: %s.\n", default_input, default_output);
			input_name = default_input;
			output_name = default_output;
		}
	}
	else {
		input_name = argv[1];
		output_name = argv[2];
	}

	start_timer(&total);
	start_timer(&stats);

	// Load file statistics and find statistics for the training data
	status = load_file(input_name, &training_file, 1000000);
	if (status != LF_ERROR_NONE) {
		printf("load_file returned error: %d\n", status);
		return 1;
	}

	training_stats = alloc_conditional_pmf_list(alphabet, training_file.columns);
	qlist = generate_codebooks(&training_file, training_stats, dist, 0.5, NULL);
	columns = qlist->columns;
	stop_timer(&stats);
	start_timer(&encoding);

	printf("Stats and codebook generation took %.4f seconds\n", get_timer_interval(&stats));

	// Allocate space for the stuff we're going to use to store what is being written to a file
	line = (char *) calloc(4096, sizeof(char));
	outline = (uint8_t *) calloc(columns, sizeof(uint8_t));
	stateline = (uint8_t *) calloc(columns, sizeof(uint8_t));
	precodeline = (uint8_t *) calloc(columns, sizeof(uint8_t));
	bits = (uint8_t *) calloc(columns, sizeof(uint8_t));
	
	// Now, open the input file and the output file
	fp = fopen(input_name, "rt");
	fbitout = fopen(output_name, "wb");
//	fref = fopen("ref.txt", "wt");
	if (!fp) {
		perror("Unable to open input file");
		exit(1);
	}
	if (!fbitout) {
		perror("Unable to open output file");
		exit(1);
	}

	// Initialize WELL state vector with libc rand (this initial vector needs to be copied to the decoder)
	srand(time(0));
	for (s = 0; s < 32; ++s) {
		qlist->well.state[s] = rand();
		// Testing with fixed state to look for consistency!
		//cb_list.well.state[s] = 0x55555555;
	}
	
	// Write the initial WELL state vector to the file first (fixed size of 32 bytes)
	fwrite(qlist->well.state, sizeof(uint32_t), 32, fbitout);

	bits_used = 0;
	j = 0;
	fgets(line, columns+2, fp);
	do {
		// State encoding is offset from the main loop iteration by one because we need
		// to use the un-encoded version to select the next codebook. After adding the precodeline
		// variable this could be fixed but now I'm too lazy to go back. Also, in a fast implementation
		// that didn't need to immediately decode to check its work, there is no point for the precode (comparison)
		// variable, so we'd have to stick to this type of loop organization

		// Select first column's codebook with no left context
		q = choose_quantizer(qlist, 0, 0);

		// Quantize and calculate error simultaneously
		// Note that in this version the quantizer outputs are 0-41, so the +33 offset is different from before
		outline[0] = q->q[line[0]-33];
		precodeline[0] = outline[0];
		error = (line[0] - outline[0] - 33)*(line[0] - outline[0] - 33);
		state_enc = find_state_encoding(q, outline[0]);
		bits[0] = cb_log2(q->output_alphabet->size);

		for (s = 1; s < columns; ++s) {
			// Quantize and compute error for MSE
			q = choose_quantizer(qlist, s, outline[s-1]);
			outline[s] = q->q[line[s]-33];
			precodeline[s] = outline[s]+33;
			error += (line[s] - outline[s] - 33)*(line[s] - outline[s] - 33);

			// Save the state encoded version over the previous value
			outline[s-1] = state_enc;
			// Bits will not be used once the arithmetic coder is added so we can afford to simply calculate them every time here for now
			bits[s] = cb_log2(q->output_alphabet->size);
			state_enc = find_state_encoding(q, outline[s]);
		}
		outline[columns-1] = state_enc;
		distortion += error / ((double) columns);

		// Bit-pack the state encoding output
		memset(stateline, 0, columns*sizeof(char));
		u = 0;
		k = 0;
		for (s = 0; s < columns; ++s) {
			// k represents the bit offset we have here
			mask = (1 << bits[s]) - 1;
			stateline[u] |= ((outline[s] & mask) << k);
			k += bits[s];
			
			// Check if there are leftover bits, put them into the next symbol
			if (k > 7) {
				u += 1;
				k -= 8;
				mask = (1 << k) - 1;
				stateline[u] |= (outline[s] >> (bits[s] - k)) & mask;
			}

			bits_used += bits[s];
		}

		// Write to file
		if (k > 0)
			u += 1;
		fwrite(stateline, sizeof(char), u, fbitout);
//		fwrite(outline, sizeof(char), columns, fbitout);
//		fwrite(precodeline, sizeof(char), columns, fref);
//		fwrite("\n", 1, 1, fref);

		// Get next line from file
		j += 1;
		fgets(line, columns+2, fp);
	} while (!feof(fp));
	fclose(fp);
	fclose(fbitout);

	distortion = distortion / ((double) j);
	stop_timer(&encoding);
	stop_timer(&total);

	// Compress with LZMA and obtain file information to calculate size and rate
	printf("MSE: %f\n", distortion);
//	printf("Filesize: %ld bytes\n", finfo.st_size);
	printf("Lines: %d\n", j);
	printf("Total bits used: %d\n", bits_used);
//	printf("Rate: %f\n", finfo.st_size*8./(((double)j)*columns));
	printf("Encoding took %.4f seconds.\n", get_timer_interval(&total));
	printf("Total time elapsed: %.4f seconds.\n", get_timer_interval(&total));

#ifndef LINUX
	system("pause");
#endif

	return 0;
}
