#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "pmf.h"
#include "util.h"

/**
 * Allocates the memory for an alphabet structure and fills the symbols
 * with a default list of 0 through size-1
 */
struct alphabet_t *alloc_alphabet(uint32_t size) {
	symbol_t i;
	struct alphabet_t *rtn = (struct alphabet_t *) calloc(1, sizeof(struct alphabet_t));
	rtn->size = size;
	rtn->symbols = (symbol_t *) calloc(size, sizeof(symbol_t));

	for (i = 0; i < size; ++i) {
		rtn->symbols[i] = i;
	}
	return rtn;
}

/**
 * Allocates a PMF structure for the given alphabet, but it does not copy the alphabet
 */
struct pmf_t *alloc_pmf(const struct alphabet_t *alphabet) {
	struct pmf_t *rtn = (struct pmf_t *) calloc(1, sizeof(struct pmf_t));
	rtn->alphabet = alphabet;
	rtn->pmf = (double *) calloc(alphabet->size, sizeof(double));
	rtn->counts = (uint32_t *) calloc(alphabet->size, sizeof(uint32_t));
	return rtn;
}

/**
 * Frees an alphabet
 */
void free_alphabet(struct alphabet_t *alphabet) {
	free(alphabet->symbols);
	free(alphabet);
}

/**
 * Frees a PMF
 */
void free_pmf(struct pmf_t *pmf) {
	free(pmf->counts);
	free(pmf);
}

/**
 * Gets the probability for a specific location, triggering lazy re-eval if
 * necessary
 */
double get_probability(struct pmf_t *pmf, uint32_t idx) {
	if (!pmf->pmf_ready)
		recalculate_pmf(pmf);
	return pmf->pmf[idx];
}

/**
 * Gets the probability for a specific symbol, triggering lazy re-eval if
 * necessary
 */
double get_symbol_probability(struct pmf_t *pmf, symbol_t symbol) {
	uint32_t idx = get_symbol_index(pmf->alphabet, symbol);

	if (!pmf->pmf_ready)
		recalculate_pmf(pmf);
	if (idx != ALPHABET_SYMBOL_NOT_FOUND)
		return pmf->pmf[idx];
	return 0.0;
}

/**
 * Calculate the entropy of this pmf in bits
 */
double get_entropy(struct pmf_t *pmf) {
	double entropy = 0.0;
	uint32_t i = 0;

	if (!pmf->pmf_ready)
		recalculate_pmf(pmf);

	for (i = 0; i < pmf->alphabet->size; ++i) {
		if (pmf->pmf[i] > 0.0) {
			entropy -= pmf->pmf[i] * log2(pmf->pmf[i]);
		}
	}

	return entropy;
}

/**
 * Calculates the Kullbeck-Leibler Divergence between two PMFs, p and q, as D(p||q)
 */
double get_kl_divergence(struct pmf_t *p, struct pmf_t *q) {
	double d = 0.0;
	uint32_t i;
	
	if (p->alphabet != q->alphabet)
		return NAN;

	if (!p->pmf_ready)
		recalculate_pmf(p);
	if (!q->pmf_ready)
		recalculate_pmf(q);

	for (i = 0; i < p->alphabet->size; ++i) {
		if (q->pmf[i] > 0) {
			if (p->pmf[i] > 0) {
				d += p->pmf[i] * log2(p->pmf[i] / q->pmf[i]);
			}
		}
	}

	return d;
}

/**
 * Combine two PMFs with two weight parameters to scale each before adding. This
 * operates based on the probabilities, not the counts, so it is suitable for use
 * in calculating the law of total probability: p(a)p(X|Y=a) + p(b)p(X|Y=b) when
 * the empirical distributions do not contain the same number of observations
 */
struct pmf_t *combine_pmfs(struct pmf_t *a, struct pmf_t *b, double weight_a, double weight_b, struct pmf_t *result) {
	uint32_t i;

	if (a->alphabet != b->alphabet || a->alphabet != result->alphabet)
		return NULL;
	
	if (!a->pmf_ready)
		recalculate_pmf(a);
	if (!b->pmf_ready)
		recalculate_pmf(b);

	for (i = 0; i < a->alphabet->size; ++i) {
		result->pmf[i] = weight_a * a->pmf[i] + weight_b * b->pmf[i];
	}
	result->pmf_ready = 1;
	return result;
}

/**
 * When counting symbols, this handles incrementing everything for the given
 * index
 */
void pmf_increment(struct pmf_t *pmf, uint32_t index) {
	pmf->counts[index] += 1;
	pmf->total += 1;
}

/**
 * Recalculates the PMF as a series of doubles from the empirical counts and total
 */
void recalculate_pmf(struct pmf_t *pmf) {
	uint32_t i;
	for (i = 0; i < pmf->alphabet->size; ++i) {
		pmf->pmf[i] = ((double) pmf->counts[i]) / pmf->total;
	}
	pmf->pmf_ready = 1;
}

/**
 * Converts a PMF that is stored as a series of doubles back to the counts representation,
 * or alternatively this can be viewed as quantizing it into a fixed point representation in
 * 0.m format
 */
void pmf_to_counts(struct pmf_t *pmf, uint32_t m) {
	uint32_t i;
	double scale = ((1 << m) - 1);

	pmf->total = 0;
	for (i = 0; i < pmf->alphabet->size; ++i) {
		pmf->counts[i] = (uint32_t) (pmf->counts[i] * scale);
		pmf->total += pmf->counts[i];
	}
}

/**
 * Looks up the index of a symbol in the given alphabet, which may be useful
 * if the alphabet doesn't start at zero, has gaps, etc.
 */
uint32_t get_symbol_index(const struct alphabet_t *alphabet, symbol_t symbol) {
	uint32_t idx;
	for (idx = 0; idx < alphabet->size; ++idx) {
		if (alphabet->symbols[idx] == symbol)
			return idx;
	}
	return ALPHABET_SYMBOL_NOT_FOUND;
}

/**
 * Displays an alphabet as "(index): 'character' <number>" one per line
 */
void print_alphabet(const struct alphabet_t *alphabet) {
	uint32_t i;
	for (i = 0; i < alphabet->size; ++i) {
		printf("(%d): '%c' <%d>\n", i, alphabet->symbols[i], alphabet->symbols[i]);
	}
}

/**
 * Displays a PMF
 */
void print_pmf(struct pmf_t *pmf) {
	uint32_t i;

	if (!pmf->pmf_ready)
		recalculate_pmf(pmf);

	for (i = 0; i < pmf->alphabet->size; ++i) {
		printf("<%d>: %.5f (%d/%d)\n", pmf->alphabet->symbols[i], pmf->pmf[i], pmf->counts[i], pmf->total);
	}
}
