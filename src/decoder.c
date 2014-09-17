
#include "util.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <time.h>

#include "codebook.h"

/**
 * Decodes the file given to produce a quality value file
 */
int main_decoder(int argc, char **argv) {
	FILE *fin, *fout;
	uint32_t columns;
	uint8_t *data, *line;
	uint8_t done = 0, next_done = 0;
	uint32_t u, k, s, j;
	uint32_t threshold = 2048;
	uint32_t bits;
	uint8_t mask;
	struct hrtimer_t timer;
	struct cond_quantizer_list_t *qlist;
	struct quantizer_t *q;
	struct alphabet_t *A = alloc_alphabet(41);

	// Check arguments
	if (argc != 4) {
		printf("Decodes a file back into a series of quality scores.\n");
		printf("Usage: %s [codebook] [input file] [output file]\n", argv[0]);
		exit(1);
	}

	start_timer(&timer);

	// Read in the codebook and find out how many columns we have
	qlist = read_codebook(argv[1], A);
	columns = qlist->columns;

	// Allocate space for decoding buffers
	data = (uint8_t *) calloc(4096, sizeof(uint8_t));
	line = (uint8_t *) calloc(columns+1, sizeof(uint8_t));
	line[columns] = '\n';

	// Open the files we're reading/writing and bail on a failure
	fin = fopen(argv[2], "rb");
	fout = fopen(argv[3], "wt");
	if (!fin) {
		perror("Unable to open input file");
		exit(1);
	}
	if (!fout) {
		perror("Unable to open output file");
		exit(1);
	}

	// First, read in the WELL state and set up the PRNG
	fread(qlist->well.state, sizeof(uint32_t), 32, fin);

	// Now, read bytes from the file and do the decoding
	fread(data, sizeof(uint8_t), 4096, fin);
	u = 0;
	j = 0;
	while (!done) {
		// Data is broken into two halves. The first 2048 bytes should be decoded, and the second half
		// is for spillover to finish whatever line we're on, before it gets moved to the front half
		// and a new spillover chunk is read in from the file

		// Handle the first column in a line specially
		q = choose_quantizer(qlist, 0, 0);
		bits = cb_log2(q->output_alphabet->size);
		mask = (1 << bits) - 1;
		line[0] = q->output_alphabet->symbols[data[u] & mask] + 33;
		k = bits;

		// If this symbol was a whole byte, advance a byte in our data stream
		if (bits == 8) {
			u += 1;
			k = 0;
		}

		// Now iterate over the rest of the columns and decode based on history
		for (s = 1; s < columns; ++s) {
//			if (s == 2)
//				printf(".");

			q = choose_quantizer(qlist, s, line[s-1]-33);
			bits = cb_log2(q->output_alphabet->size);

			// Unpack bits from the current byte
			mask = (1 << bits) - 1;
			line[s] = ((((int) data[u]) & 0xff) >> k) & mask;
			k += bits;

			// If there are stray bits we need from the next byte, grab those too
			if (k > 8) {
				u += 1;
				mask = (1 << (k - 8)) - 1;
				line[s] |= (data[u] & mask) << (8 - k + bits);
				k -= 8;
			}

			// Undo state encoding
			line[s] = q->output_alphabet->symbols[line[s]] + 33;
		}

		// Write this line to the output file, note \n at the end of the line buffer to get the right length
		fwrite(line, columns+1, sizeof(uint8_t), fout);
		j += 1;

		// We always move to the next byte after a complete line
		u += 1;

		// Check if we need to shift the buffers and read more data
		if (u >= threshold) {
			if (!next_done) {
				memmove(data, &data[u], 4096 - u);
				threshold = fread(&data[4096 - u], sizeof(uint8_t), u, fin);
				// EOF is only triggered when we can't pull in anything at all
				if (feof(fin)) {
					next_done = 1;
					threshold = 4096 - u + threshold - 1;
				}
				else {
					threshold = 2048;
				}
				// Start over near the beginning of the buffer
				u = 0;
			}
			else {
				done = 1;
			}
		}
	}

	// Close the files and we are done decoding
	fclose(fin);
	fclose(fout);

	stop_timer(&timer);
	printf("Decoded %d lines in %f seconds.\n", j, get_timer_interval(&timer));

#ifndef LINUX
	system("pause");
#endif

	return 0;
}
