#include "codebook.h"

// Legacy stuff for the current version; this should be updated as things are implemented
// using the new C-based calculations

/**
 * Reads in the codebook in the specified filename, calculates how many columns it is configured
 * for, and prepares the codebook list structure necessary to use for encoding with it
 */
uint32_t read_codebook(const char *filename, struct codebook_list_t *cb_list, uint8_t symbols) {
	FILE *fp;
	uint32_t columns;
	uint32_t column, j;
	char line[MAX_CODEBOOK_LINE_LENGTH];

	fp = fopen(filename, "rt");
	if (!fp) {
		perror("Unable to open codebook file");
		exit(1);
	}

	// Figure out how many columns the data has, accounting for the newline at the end of the line
	fgets(line, MAX_CODEBOOK_LINE_LENGTH, fp);
	columns = strlen(line) - 1;

	// Initialize codebook based on knowing how many lines we have
	init_codebook_list(cb_list, symbols, columns);

	// Skip the first line (already read) and second line because we don't actually care about how many states there are
	fgets(line, MAX_CODEBOOK_LINE_LENGTH, fp);

	// Next line is the selection ratio between the two codebooks
	fgets(line, MAX_CODEBOOK_LINE_LENGTH, fp);
	for (j = 0; j < columns; ++j) {
		cb_list->ratio[j] = line[j] - 33;
	}

	// Now, the lines in file alternate definitions for each codebook (low then high), so copy them over and generate the corresponding unique lists
	fgets(line, MAX_CODEBOOK_LINE_LENGTH, fp);
	memcpy(cb_list->low[0][0].quantizer, line, symbols * sizeof(uint8_t));
	generate_uniques(&cb_list->low[0][0]);

	fgets(line, MAX_CODEBOOK_LINE_LENGTH, fp);
	memcpy(cb_list->high[0][0].quantizer, line, symbols * sizeof(uint8_t));
	generate_uniques(&cb_list->high[0][0]);

	// Parse remaining lines as code books conditional on previous column values
	for (column = 1; column < columns; ++column) {
		fgets(line, MAX_CODEBOOK_LINE_LENGTH, fp);
		for (j = 0; j < symbols; ++j) {
			memcpy(cb_list->low[column][j].quantizer, &line[j*symbols], symbols * sizeof(uint8_t));
			generate_uniques(&cb_list->low[column][j]);
		}

		fgets(line, MAX_CODEBOOK_LINE_LENGTH, fp);
		for (j = 0; j < symbols; ++j) {
			memcpy(cb_list->high[column][j].quantizer, &line[j*symbols], symbols * sizeof(uint8_t));
			generate_uniques(&cb_list->high[column][j]);
		}
	}

	// Done reading code book
	fclose(fp);
	return columns;
}

/**
 * Initializes a codebook list with storage space for an adequate number of codebooks
 * all initialized to all zeros, for the given number of possible symbols
 */
void init_codebook_list(struct codebook_list_t *list, uint8_t symbols, uint32_t columns) {
	
	// First, allocate arrays based on number of columns within the codebook list
	list->high = (struct codebook_t **) calloc(columns, sizeof(struct codebook_t *));
	list->low = (struct codebook_t **) calloc(columns, sizeof(struct codebook_t *));
	list->ratio = (uint8_t *) calloc(columns, sizeof(uint8_t));
	list->select_count = (uint32_t *) calloc(columns, sizeof(uint32_t));
	list->symbols = symbols;
	list->columns = columns;

	// Now allocate for each array
	init_codebook_array(list->high, symbols, columns);
	init_codebook_array(list->low, symbols, columns);

	// Make sure that WELL starts off reasonable
	list->well.n = 0;
}

/**
 * Initializes an array of codebooks, used when initializing a codebook list
 */
void init_codebook_array(struct codebook_t **cb, uint8_t symbols, uint32_t columns) {
	uint32_t c, s;

	// First column is special in that it only has one codebook, because there is no left context
	cb[0] = (struct codebook_t *) calloc(1, sizeof(struct codebook_t));
	cb[0]->quantizer = (uint8_t *) calloc(symbols, sizeof(uint8_t));
	cb[0]->uniques = (uint8_t *) calloc(symbols, sizeof(uint8_t));
	cb[0]->symbols = symbols;

	for (c = 1; c < columns; ++c) {
		cb[c] = (struct codebook_t *) calloc(symbols, sizeof(struct codebook_t));
		for (s = 0; s < symbols; ++s) {
			cb[c][s].symbols = symbols;
			cb[c][s].quantizer = (uint8_t *) calloc(symbols, sizeof(uint8_t));
			cb[c][s].uniques = (uint8_t *) calloc(symbols, sizeof(uint8_t));
		}
	}
}

/**
 * Walks over the quantizer string to determine how many unique symbols are present
 */
void generate_uniques(struct codebook_t *cb) {
	uint8_t u = 0;
	uint8_t s;

	cb->uniques[0] = cb->quantizer[0];
	for (s = 1; s < cb->symbols; ++s) {
		if (cb->quantizer[s] != cb->uniques[u]) {
			u += 1;
			cb->uniques[u] = cb->quantizer[s];
		}
	}

	cb->actual_unique_count = u+1;
	cb->bits = cb_log2(cb->actual_unique_count);
}

/**
 * Selects a codebook for the given column from the codebook list with the appropriate ratio
 */
struct codebook_t *choose_codebook(struct codebook_list_t *list, uint32_t column, uint8_t prev_value) {
	if (well_1024a(&list->well) % 100 >= list->ratio[column]) {
		list->select_count[column] += 1;
		return &list->high[column][prev_value];
	}
	return &list->low[column][prev_value];
}

/**
 * Converts the quality score into a state number that can be stored in fewer bits
 * @param value is an ascii character to be converted into a numbered state
 */
uint8_t find_state_encoding(struct codebook_t *codebook, uint8_t value) {
	uint8_t u;

	for (u = 0; u < codebook->actual_unique_count; ++u) {
		if (codebook->uniques[u] == value)
			return u;
	}
	return u;
}

/**
 * Displays a codebook on STDOUT
 */
void print_codebook(struct codebook_t *cb) {
	uint8_t s = 0;
	uint8_t *tmp;

	tmp = (uint8_t *) calloc(cb->symbols+1, sizeof(uint8_t));

	for (s = 0; s < cb->symbols; ++s) {
		memcpy(tmp, cb[s].quantizer, cb->symbols);
		printf("%d (%c):\t%s\n", s, (uint8_t) (s+33), tmp);
	}

	free(tmp);
}

