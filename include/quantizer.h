#ifndef _QUANTIZER_H_
#define _QUANTIZER_H_

#include <stdint.h>

#include "pmf.h"
#include "distortion.h"
#include "util.h"

#define QUANTIZER_MAX_ITER		100

/**
 * Structure holding information about a quantizer, which just maps input symbols
 * to output symbols for a specific alphabet
 */
struct quantizer_t {
	const struct alphabet_t *alphabet;
	struct alphabet_t *output_alphabet;
	uint8_t *q;
};

// Memory management
struct quantizer_t *alloc_quantizer(const struct alphabet_t *);
void free_quantizer(struct quantizer_t *);

// Generates a quantizer via optimization
struct quantizer_t *generate_quantizer(struct pmf_t *pmf, struct distortion_t *dist, uint8_t states, double *dist_out);

// Calculate the output pmf when the quantizer is applied to the input pmf
void apply_quantizer(struct quantizer_t *q, struct pmf_t *restrict pmf, struct pmf_t *restrict output);

// Display/debugging
void print_quantizer(struct quantizer_t *);

#endif
