#include "codebook.h"
#include "pmf.h"
#include "lines.h"

/**
 * To compute stats for the training data, we will need a set of conditional PMFs, one
 * per column
 * @param alphabet The symbol alphabet for each column
 * @param columns The number of columns to allocate conditional PMFs for
 */
struct cond_pmf_list_t *alloc_conditional_pmf_list(const struct alphabet_t *alphabet, uint32_t columns) {
	uint32_t count = 1 + alphabet->size*(columns-1);
	struct cond_pmf_list_t *list = (struct cond_pmf_list_t *) calloc(1, sizeof(struct cond_pmf_list_t));

	// We need one array of PMF pointers that will index into the buffer allocated above, for the columns
	list->columns = columns;
	list->alphabet = alphabet;
	list->pmfs = (struct pmf_t **) calloc(count, sizeof(struct pmf_t *));

	// All PMFs are stored in a flat array, the accessor function will resolve a PMF's address
	for (i = 0; i < count; ++i) {
		list->pmfs[i] = alloc_pmf(alphabet);
	}

	return list;
}

/**
 * Deallocate the PMF list given and unallocate the two allocated memory blocks
 * @param list The conditional pmf list to deallocate
 */
void free_conditional_pmf_list(struct cond_pmf_list_t *list) {
	uint32_t count = 1 + list->alphabet->size * (list->columns - 1);
	uint32_t i;

	for (i = 0; i < count; ++i) {
		free_pmf(list->pmfs[0]);
	}
	free(list);
}

/**
 * Allocate the quantizer list structure and the first level of array based on columns
 * @param columns The number of columns for which we have quantizers
 * @return Pointer to conditional quantizer list structure
 */
struct cond_quantizer_list_t *alloc_conditional_quantizer_list(uint32_t columns) {
	struct cond_quantizer_list_t *rtn = (struct cond_quantizer_list_t *) calloc(1, sizeof(struct cond_quantizer_list_t));
	rtn->columns = columns;
	rtn->input_alphabets = (struct alphabet_t **) calloc(columns, sizeof(struct alphabet_t *));
	rtn->q = (struct quantizer_t ***) calloc(columns, sizeof(struct quantizer_t **));
	return rtn;
}

/**
 * Deallocate the quantizer list as well as any alphabets or pmfs that are stored
 * @param list The conditional quantizer list to deallocate
 */
void free_cond_quantizer_list(struct cond_quantizer_list_t *list) {
	uint32_t i, j;

	for (i = 0; i < columns; ++i) {
		if (list->q[i]) {
			for (j = 0; j < list->input_alphabets[i].size; ++i) {
				if (list->q[i][j])
					free_pmf(list->q[i][j]);
			}
			free_alphabet(list->alphabets[i]);
			free(list->q[i]);
		}
	}

	free(list->q);
	free(list->input_alphabets);
	free(list);
}

/**
 * Initialize the information within a quantizer for the given column. This can't be done
 * at allocation time because we don't know everything about this column until we get here
 * during the optimization process
 * @param list The conditional quantizer list to update
 * @param column The column to initialize
 * @param input_union The alphabet of all possible left context symbols
 */
void cond_quantizer_init_column(struct cond_quantizer_list_t *list, uint32_t column, const struct alphabet_t *input_union) {
	list->input_alphabets[column] = duplicate_alphabet(input_union);
	list->q[column] = (struct quantizer_t **) calloc(input_union->size, sizeof(struct quantizer_t *));
}

/**
 * Find a PMF for a specific column with the specific previous value
 */
struct pmf_t *get_cond_pmf(struct cond_pmf_list_t *list, uint32_t column, symbol_t prev) {
	if (column == 0)
		return list->pmfs[0];
	return list->pmfs[1 + (column-1)*list->alphabet->size + prev];
}

/**
 * Get a quantizer by its indexed location within the quantizer list for a column
 */
struct quantizer_t *get_cond_quantizer_indexed(struct cond_quantizer_list_t *list, uint32_t column, uint32_t index) {
	return list->q[column][index];
}

/**
 * Get a quantizer by its left context symbol
 */
struct quantizer_t *get_cond_quantizer(struct cond_quantizer_list_t *list, uint32_t column, symbol_t prev) {
	uint32_t idx = get_symbol_index(list->input_alphabets[column], prev);
	if (idx != ALPHABET_SYMBOL_NOT_FOUND)
		return get_cond_quantizer_indexed(list, column, idx);
	return NULL;
}

/**
 * Stores the given quantizer at the appropriate index corresponding to the left context symbol given
 * for the specific column
 */
void store_cond_quantizer(struct cond_quantizer_t *q, struct cond_quantizer_list_t *list, uint32_t column, symbol_t prev) {
	uint32_t idx = get_symbol_indx(list->input_alphabets, prev);
	list->q[column][idx] = q;
}

/**
 * Given a quality info struct, which is assumed to have loaded the training set, and a
 * set of output conditional pmf structures, calculate the statistics of the data
 */
void calculate_statistics(struct quality_info_t *info, struct cond_pmf_list_t *pmf_list) {
	uint32_t block_idx, line_idx;
	struct line_t *line;

	for (block_idx = 0; block_idx < info->block_count; ++block_idx) {
		for (line = 0; line < info->blocks[block_idx].count; ++line) {
			line = &info->blocks[block_idx].lines[line_idx];
			pmf_increment(get_cond_pmf(pmf_list, 0, 0), line->data[0]);
			for (column = 1; column < info->columns; ++column) {
				pmf_increment(get_cond_pmf(pmf_list, column, line->data[column-1]), line->data[column]);
			}
		}
	}
}

/**
 * Calculates the integer number of states to use for each column according to the estimate
 * of conditional entropy from the baseline statistics
 */
void find_bit_allocation(struct cond_pmf_list_t *pmf_list, double comp, double **high, double **low, double **ratio, uint32_t mode) {
	*high = (double *) calloc(pmf_list->columns, sizeof(double));
	*low = (double *) calloc(pmf_list->columns, sizeof(double));
	*ratio = (double *) calloc(pmf_list->columns, sizeof(double));
	double *entropies = (double *) _alloca(pmf_list->columns*sizeof(double));
	double h_lo, h_hi;
	uint32_t i, j;
	struct pmf_list_t *uc_pmf_list = alloc_pmf_list(pmf_list->columns, pmf_list->alphabet);
	struct pmf_t *pmf_temp;

	// Column 0 pmf is the unconditional one, copy it to the unconditional pmf list
	pmf_temp = get_cond_pmf(pmf_list, 0, 0);
	combine_pmfs(pmf_temp, pmf_temp, 1.0, 0.0, &uc_pmf_list->pmfs[0]);

	// Find unconditional pmfs for each column remaining
	for (i = 1; i < pmf_list->columns; ++i) {
		for (j = 0; j < pmf_list->alphabet->size; ++j) {
			pmf_temp = get_cond_pmf(pmf_list, i, j);
			combine_pmfs(&uc_pmf_list->pmfs[i], pmf_temp, 1.0, get_probability(&uc_pmf_list->pmfs[i-1], j), &uc_pmf_list->pmfs[i]);
		}
	}

	// Alloca doesn't zero memory for us
	memset(entropies, 0, pmf_list->columns*sizeof(double));

	// Column 0 is handled specially because it only has one left context
	entropies[0] = get_entropy(get_cond_pmf(pmf_list, 0, 0));

	// Rest of the columns are handled identically
	for (i = 1; i < pmf_list->columns; ++i) {
		pmf_temp = &uc_pmf_list->pmfs[i-1];
		for (j = 0; j < pmf_list->alphabet->size; ++j) {
			entropies[i] += get_probability(pmf_temp, j) * get_entropy(get_cond_pmf(pmf_list, i, j));
		}
	}

	// Compute number of states used based on mode parameter
	// r = ratio:
	// H = rH_lo + (1-r)H_hi
	// H - H_hi = r(H_lo - H_hi)
	// r = (H - H_hi) / (H_lo - H_hi)
	for (i = 0; i < pmf_list->columns; ++i) {
		switch(mode) {
			case BIT_ALLOC_MODE_INT_STATES:
				h_lo = pow(2, entropies[i]);
				low[i] = floor(h_lo);
				high[i] = ceil(h_lo);
				h_lo = log2(low[i]);
				h_hi = log2(high[i]);
				ratio[i] = (entropies[i] - h_hi) / (h_lo - h_hi);
				break;
			case BIT_ALLOC_MODE_INT_POWER:
				h_lo = floor(entropies[i]);
				h_hi = ceil(entropies[i]);
				low[i] = pow(2, h_lo);
				high[i] = pow(2, h_hi);
				ratio[i] = (entropies[i] - h_hi) / (h_lo - h_hi);
				break;
			case BIT_ALLOC_MODE_NO_MIX:
			default:
				ratio[i] = 1;
				low[i] = floor(pow(2, entropies[i]));
				high[i] = 0.0;
				break;
		}
	}

	free_pmf_list(uc_pmf_list);
}

/**
 * Given the statistics calculated before, we need to compute the entire codebook's worth of
 * quantizers, as well as all of the PMFs and related stats
 * TODO: Add hi states to generation
 */
struct cond_quantizer_list_t *generate_codebooks(struct quality_info_t *info, struct cond_pmf_list_t *in_pmfs, struct distortion_t *dist, double comp, uint32_t mode) {
	// Stuff for state allocation and mixing
	double *hi_states, *lo_states, *state_ratio;

	// Miscellaneous variables
	uint32_t column, i, j, k;
	symbol_t x, q;
	double weight, norm;

	// Output list of conditional quantizers
	struct cond_quantizer_list_t *q_list;

	// Constant alphabet of all possible input symbols
	const struct alphabet_t *A = in_pmfs->alphabet;

	// Temporary/extra pointers
	struct quantizer_t *q_temp;
	struct pmf_t *pmf_temp;

	// List of conditionally quantized PMFs, conditioned on which quantizer was chosen
	struct pmf_list_t **cq_pmf_list;
	
	// List of conditionally quantized PMFs after quantizer has been added out
	struct pmf_list_t *xpmf_list;

	// List of conditionally quantized PMFs after the next quantizer was applied
	struct pmf_list_t *qpmf_list;

	// "Alphabet" of output quantizers for the previous column and PMF over that set
	struct alphabet_t *q_alphabet;
	struct pmf_t *q_pmf;

	// Alphabet of all possible quantizer outputs from the previous column
	struct alphabet_t *q_output_union;

	// First we need to know what our training stats are and figure out how many states
	// to put into each quantizer
	calculate_statistics(info, in_pmfs);
	find_bit_allocation(pmf_list, comp, &hi_states, &lo_states, &state_ratio, mode);

	// For the first column, quantizer isn't conditional, so find it directly
	q_temp = generate_quantizer(get_cond_pmf(in_pmfs, 0, 0), dist, lo_states[0], NULL);
	store_cond_quantizer(q_temp, q_list, 0, 0);

	// Set up a quantizer alphabet for column zero. The cond quantizer list duplicates the
	// alphabet internally so we don't need to worry about duplicating it
	q_alphabet = alloc_alphabet(1);
	cond_quantizer_init_column(q_list, 0, q_alphabet);

	// Initialize a 100% PMF for this quantizer alphabet
	q_pmf = alloc_pmf(q_alphabet);
	pmf_increment(q_pmf, 0);
	q_output_union = duplicate_alphabet(q_temp->output_alphabet);

	// Previous qpmf list needs to be initialized for a single quantizer's output
	prev_qpmf_list = alloc_pmf_list(1, A);
	apply_quantizer(q_temp, get_cond_pmf(in_pmfs, 0, 0), &prev_qpmf_list->pmfs[0]);

	// Iterate over remaining columns and compute the quantizers and quantized PMFs
	for (column = 1; column < info->columns; ++column) {
		// Per quantizer from previous column, we'll have a list of conditional quantized PMFs
		cq_pmf_list = (struct pmf_list_t **) calloc(q_alphabet->size, sizeof(struct pmf_list_t *));
		for (i = 0; i < q_alphabet->size; ++i) {
			// Allocate a complete set of PMFs, to simplify lookup later
			q_temp = get_cond_quantizer_indexed(q_list, column-1, i);
			cq_pmf_list[i] = alloc_pmf_list(A->size, A);

			// For each conditional symbol in the unquantized input, find corresponding quantized value
			// and then add probability mass to that pmf
			for (j = 0; j < A->size; ++j) {
				q = q_temp->q[j];
				pmf_temp = &cq_pmf_list[i]->pmfs[q];
				combine_pmfs(pmf_temp, get_cond_pmf(in_pmfs, column, j), 1.0, 1.0, pmf_temp);
			}
		}

		// Now, based on the probability of each quantizer, we need to merge the conditional
		// quantized pmfs given quantizer into a single list of conditional pmfs depending on
		// only the terms in the output union alphabet
		xpmf_list = alloc_pmf_list(A->size, A);
		for (j = 0; j < q_output_union->size; ++j) {
			q = q_output_union->symbols[j];
			norm = 0;

			// First, find the normalizing constant by adding probabilities where this symbol
			// is produced
			for (i = 0; i < q_alphabet->size; ++i) {
				q_temp = get_cond_quantizer(q_list, column-1, i);
				if (alphabet_contains(q_temp->output_alphabet, q)) {
					norm += get_probability(q_pmf, i);
				}
			}

			// Then, go over the quantizer alphabet again and merge PMFs weighted by the probability
			// of each quantizer and normalized by the total probability of all quantizers producing
			// a given symbol
			for (i = 0; i < q_alphabet->size; ++i) {
				q_temp = get_cond_quantizer(q_list, column-1, i);
				if (alphabet_contains(q_temp->output_alphabet, q)) {
					weight = get_probability(q_pmf, i);
					pmf_temp = &xpmf_list->pmfs[j];
					combine_pmfs(pmf_temp, &cq_pmf_list[i]->pmfs[q], 1.0, weight/norm, pmf_temp);
				}
			}
		}

		// So, xpmf_list is now the set of all PMFs for this column conditioned on the previous
		// column's quantized value, with history, choice of quantizer, etc. all added out.
		// Now we need to compute the stats that we need for the next iteration

		// So, first, compute the next set of quantizers and associated output PMFs
		cond_quantizer_init_column(q_list, column, q_output_union);
		qpmf_list = alloc_pmf_list(q_output_union->size, A);
		for (j = 0; j < q_output_union->size; ++j) {
			q = q_output_union->symbols[j];

			// Find and save quantizer
			q_temp = generate_quantizer(&xpmf_list->pmfs[q], dist, lo_states[column], NULL);
			store_cond_quantizer(q_temp, q_list, column, q);

			// Find the PMF of the quantizer's output
			apply_quantizer(q_temp, &xpmf_list->pmfs[q], &qpmf_list->pmfs[j]);
		}

		// Compute the next output alphabet union over all quantizers for this column
		next_output_union = duplicate_alphabet(get_cond_quantizer_indexed(q_list, column, 0)->output_alphabet);
		for (j = 1; j < q_output_union->size; ++j) {
			alphabet_union(next_output_union, get_cond_quantizer_indexed(q_list, column, j)->output_alphabet, next_output_union);
		}

		// Find the pmf of choice of quantizer for this column
		next_q_alphabet = alloc_alphabet(q_output_union->size);
		next_q_pmf = alloc_pmf(next_q_alphabet);
		for (j = 0; j < next_q_alphabet->size; ++j) {
			// next_q_alphabet is the same size as q_output_union, but the first indexes quantizers
			// and the second is the unique set of quantizer outputs
			q = q_output_union->symbols[j];
			for (i = 0; i < q_alphabet->size; ++i) {
				// Total probability over PMF of each symbol conditioned on the quantizer that produced it
				// Output symbol selects which quantizer is used from this column
				next_q_pmf->pmf[j] += get_probability(q_pmf, i) * get_probability(&prev_qpmf_list->pmfs[i], q);
			}
		}

		// Finally, deallocate memory that was used on this iteration
		free(q_output_union);
		q_output_union = next_output_union;

		free_pmf_list(prev_qpmf_list);
		prev_qpmf_list = qpmf_list;
		free_pmf_list(xpmf_list);
		for (i = 0; i < q_alphabet->size; ++i) {
			free_pmf_list(cq_pmf_list[i]);
		}
		free(cq_pmf_list);

		free_pmf(q_pmf);
		q_pmf = next_q_pmf;
		free_alphabet(q_alphabet);
		q_alphabet = next_q_alphabet;
	}

	// Final cleanup, things we saved at the end of the final iteration that aren't needed
	free(q_pmf);
	free_pmf_list(qpmf_list);
	free(q_alphabet);

	// TODO: Generate a codebook-format organization of the quantizers for quick lookup during encoding
	return q_list;
}

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

