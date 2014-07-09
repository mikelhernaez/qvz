#ifndef _CODEBOOK_H_
#define _CODEBOOK_H_
/**
 * Functions and definitions relating to reading codebooks from files, used
 * for both the encoder and decoder code
 */

#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "well.h"
#include "pmf.h"
#include "distortion.h"
#include "util.h"

/**
 * Stores an array of conditional PMFs for the current column given the previous
 * column. PMF pointers are stored in a flat array so don't try to find the PMF you
 * want directly--use the accessor
 */
struct cond_pmf_list_t {
	uint32_t columns;
	const struct alphabet_t *alphabet;
	struct pmf_t **pmfs;
};

/**
 * For all columns, store the relevant PMFs (regular or conditional), the quantizers,
 * the PMFs describing quantizer outputs
 */

// Memory management
struct cond_pmf_list_t *alloc_conditional_pmf_list(const struct alphabet_t *alphabet, uint32_t columns);
void free_cond_pmf_list(struct cond_pmf_list_t *);

// Accessors
struct pmf_t *get_cond_pmf(struct cond_pmf_list_t *list, uint32_t column, symbol_t prev);

// Legacy stuff to be converted still

#define MAX_CODEBOOK_LINE_LENGTH 4096

struct codebook_t {
	uint8_t *quantizer;				// Array of quantizer mapping (index is input)
	uint8_t *uniques;				// Array of unique values in quantizer list
	uint8_t max_unique_count;		// Maximum number of uniques allowed. Unused?
	uint8_t actual_unique_count;	// Actual number of unique elements
	uint8_t bits;					// Number of bits used for state encoding this codebook
	uint8_t symbols;				// Possible number of symbols in alphabet (length of quantizer)
};

struct codebook_list_t {
	struct codebook_t **high;
	struct codebook_t **low;
	uint8_t *ratio;
	uint32_t *select_count;
	uint8_t symbols;
	uint32_t columns;
	struct well_state_t well;
};

// Master function to read a codebook from a file
uint32_t read_codebook(const char *filename, struct codebook_list_t *list, uint8_t symbols);

// Initialization and parsing
void init_codebook_list(struct codebook_list_t *list, uint8_t symbols, uint32_t columns);
void init_codebook_array(struct codebook_t **cb, uint8_t symbols, uint32_t columns);
void generate_uniques(struct codebook_t *cb);

// Operational functions
struct codebook_t *choose_codebook(struct codebook_list_t *list, uint32_t column, uint8_t prev_value);
uint8_t find_state_encoding(struct codebook_t *codebook, uint8_t value);

// Debugging
void print_codebook(struct codebook_t *cb);

#endif
