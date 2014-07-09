
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
const char *default_codebook = "codebook.txt";
const char *default_input = "quality_lines.txt";
const char *default_output = "c_bitpacked_quality_lossy.txt";

/**
 * Actually does the encoding (and also tests decoding to ensure correctness)
 */
int main(int argc, char **argv) {
	FILE *fp, *fbitout; //, *fref;
	uint32_t columns;
	uint32_t bits_used;
	uint32_t j = 0, s = 0, u = 0, k = 0;
	double distortion = 0.0;
	uint32_t error = 0;
	uint8_t state_enc = 0;
	uint8_t mask = 0;
	struct codebook_list_t cb_list;
	struct codebook_t *prev_cb;
	struct _stat finfo;
	char *line, *cmd;
	uint8_t *outline, *stateline, *precodeline;
	uint8_t *bits;
	struct hrtimer_t timer;
	const char *codebook_name, *input_name, *output_name;

	// Check for arguments and get file names
	if (argc != 4) {
		// Print usage guidelines
		printf("Usage: %s [codebook] [input file] [output file]\n", argv[0]);
		
		// Attempted to run from command line but not all arguments are present
		if (argc > 1) {
			printf("Incorrect number of arguments supplied!");
			exit(1);
		}
		else {
			printf("Using default codebook: %s\nDefault input file: %s\nDefault output file: %s.\n", default_codebook, default_input, default_output);
			codebook_name = default_codebook;
			input_name = default_input;
			output_name = default_output;
		}
	}
	else {
		codebook_name = argv[1];
		input_name = argv[2];
		output_name = argv[3];
	}

	start_timer(&timer);

	// Read the codebook file and find out how many columns it is configured for
	columns = read_codebook(codebook_name, &cb_list, SYMBOLS);

	// Allocate space for the stuff we're going to use to store what is being written to a file
	line = (char *) calloc(4096, sizeof(char));
	cmd = (char *) calloc(4096, sizeof(char));
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
		cb_list.well.state[s] = rand();
		//cb_list.well.state[s] = 0x55555555;
	}
	
	// Write the initial WELL state vector to the file first (fixed size of 32 bytes)
	fwrite(cb_list.well.state, sizeof(uint32_t), 32, fbitout);

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
		prev_cb = choose_codebook(&cb_list, 0, 0);

		// Quantize and calculate error simultaneously
		outline[0] = prev_cb->quantizer[line[0]-33];
		precodeline[0] = outline[0];
		error = (line[0] - outline[0])*(line[0] - outline[0]);
		state_enc = find_state_encoding(prev_cb, outline[0]);
		bits[0] = prev_cb->bits;

		for (s = 1; s < columns; ++s) {
			// Quantize and compute error for MSE
			prev_cb = choose_codebook(&cb_list, s, outline[s-1]-33);
			outline[s] = prev_cb->quantizer[line[s]-33];
			precodeline[s] = outline[s];
			error += (line[s] - outline[s])*(line[s] - outline[s]);

			// Check for invalid encoding
			if (outline[s] == ' ') {
				printf("Codebook produced a space for column %d on line %d. Previous symbol: '%c'\n", s, j, outline[s-1]);
			}

			// Save the state encoded version over the previous value
			outline[s-1] = state_enc;
			bits[s] = prev_cb->bits;
			state_enc = find_state_encoding(prev_cb, outline[s]);
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
	stop_timer(&timer);

#ifndef LINUX
	sprintf(cmd, "\"c:\\program files\\7-zip\\7z\" a -mx9 %s.7z %s", output_name, output_name);
#else
	sprintf(cmd, "7z a -mx9 %s.7z %s", output_name, output_name);
#endif

	system(cmd);
	sprintf(cmd, "%s.7z", output_name);
	_stat(cmd, &finfo);

	// Compress with LZMA and obtain file information to calculate size and rate
	printf("MSE: %f\n", distortion);
	printf("Filesize: %ld bytes\n", finfo.st_size);
	printf("Lines: %d\n", j);
	printf("Total bits used: %d\n", bits_used);
	printf("Rate: %f\n", finfo.st_size*8./(((double)j)*columns));
	printf("Encoding (pre-LZMA) took %f seconds.\n", get_timer_interval(&timer));

	return 0;
}
