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

// Modes for bit allocation
#define BIT_ALLOC_MODE_INT_STATES		1			// Uses floor(2^H)/ceil(2^H)
#define BIT_ALLOC_MODE_INT_POWER		2			// Uses 2^floor(H)/2^ceil(H)
#define BIT_ALLOC_MODE_NO_MIX			3			// Uses only floor(2^H)

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
 * Stores an array of quantizer pointers for the column for all possible left context
 * values. Unused ones are left as null pointers. This is also stored as a flat array
 * so the accessor must be used to look up the correct quantizer
 * The dreaded triple pointer is used to store an array of (different length) arrays
 * of pointers to quantizers
 */
struct cond_quantizer_list_t {
	uint32_t columns;
	const struct alphabet_t **input_alphabets;
	struct quantizer_t ***q;
};

// Memory management
struct cond_pmf_list_t *alloc_conditional_pmf_list(const struct alphabet_t *alphabet, uint32_t columns);
struct cond_quantizer_list_t *alloc_conditional_quantizer_list(uint32_t columns);
void free_cond_pmf_list(struct cond_pmf_list_t *);
void free_cond_quantizer_list(struct cond_quantizer_list_t *);

// Per-column initializer for conditional quantizer list
void cond_quantizer_init_column(struct cond_quantizer_list_t *list, uint32_t column, const struct alphabet_t *input_union);

// Accessors
struct pmf_t *get_cond_pmf(struct cond_pmf_list_t *list, uint32_t column, symbol_t prev);
struct quantizer_t *get_cond_quantizer_indexed(struct cond_quantizer_list_t *list, uint32_t column, uint32_t index);
struct quantizer_t *get_cond_quantizer(struct cond_quantizer_list_t *list, uint32_t column, symbol_t prev);
void store_cond_quantizer(struct cond_quantizer_t *q, struct cond_quantizer_list_t *list, uint32_t column, symbol_t prev);

// Meat of the implementation
void calculate_statistics(struct quality_info_t *, struct cond_pmf_list_t *);
void find_bit_allocation(struct cond_pmf_list_t *pmf_list, double comp, double **high, double **low, double **ratio, uint32_t mode);
struct cond_quantizer_list_t *generate_codebooks(struct quality_info_t *info, struct cond_pmf_list_t *in_pmfs, struct distortion_t *dist, double comp, uint32_t mode);

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
