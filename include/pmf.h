#ifndef _PMF_H_
#define _PMF_H_

#include <stdint.h>

// Used to indicate a symbol not found during index lookup
#define ALPHABET_SYMBOL_NOT_FOUND			UINT32_MAX

// Unfortunately this is a bit brittle so don't change it
typedef symbol_t uint8_t;

/**
 * Structure that stores information about an alphabet including
 * the number of symbols and a list of the symbols themselves
 */
struct alphabet_t {
	uint32_t size;
	symbol_t *symbols;
};

/**
 * Structure for defining and storing a single PMF in a manner that is
 * useful for computing empirical PMFs, but also allows PMFs to be manipulated
 * as a set of probabilities (not just as empirical counts)
 */
struct pmf_t {
	uint8_t pmf_ready;
	const struct alphabet_t *alphabet;
	double *pmf;
	uint32_t *counts;
	uint32_t total;
};

// Memory management functions
struct alphabet_t *alloc_alphabet(uint32_t size);
struct pmf_t *alloc_pmf(const struct alphabet_t *);
void free_alphabet(struct alphabet_t *);
void free_pmf(struct pmf_t *);

// PMF access
uint32_t get_symbol_index(const struct alphabet_t *alphabet, symbol_t symbol);
double get_probability(struct pmf_t *pmf, uint32_t idx);
double get_symbol_probability(struct pmf_t *pmf, symbol_t symbol);
double get_entropy(struct pmf_t *pmf);
double get_kl_divergence(struct pmf_t *p, struct pmf_t *q);

// PMF Manipulation
struct pmf_t *combine_pmfs(struct pmf_t *a, struct pmf_t *b, double weight_a, double weight_b, struct pmf_t *result);
void pmf_increment(struct pmf_t *pmf, uint32_t index);
void recalculate_pmf(struct pmf_t *);
void pmf_to_counts(struct pmf_t *pmf, uint32_t m);

// Alphabet search
uint32_t get_symbol_index(const struct alphabet_t *alphabet, symbol_t symbol);

// Display routines
void print_alphabet(const struct alphabet_t *);
void print_pmf(struct pmf_t *);

#endif