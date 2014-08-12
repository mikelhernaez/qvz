#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "lines.h"
#include "pmf.h"
#include "distortion.h"
#include "quantizer.h"
#include "codebook.h"

int main(int argc, char **argv) {
	struct distortion_t *dist = generate_distortion_matrix(41, DISTORTION_MSE);
	struct alphabet_t *alphabet = alloc_alphabet(41);
	struct quality_file_t training_file;
	uint32_t status;
	char *input_path;
	char *output_path;
	struct hrtimer_t timer;
	struct cond_pmf_list_t *training_stats;
	uint32_t i, j, k;
	uint32_t *hi, *lo;
	double *ratio;

	if (argc != 3) {
		printf("Usage: %s [input file] [output file]\n", argv[0]);
		return 1;
	}

	input_path = argv[1];
	output_path = argv[2];

	start_timer(&timer);

	status = load_file(input_path, &training_file, 1000000);
	if (status != LF_ERROR_NONE) {
		printf("load_file returned error: %d\n", status);
		return 1;
	}

	training_stats = alloc_conditional_pmf_list(alphabet, training_file.columns);
	calculate_statistics(&training_file, training_stats);
	find_bit_allocation(training_stats, 0.5, &hi, &lo, &ratio, BIT_ALLOC_MODE_INT_STATES);

	// Debug: verify the PMFs that were calculated
	for (i = 0; i < training_file.columns; ++i) {
		printf("Column %d: %d / %d (%f)\n", i, hi[i], lo[i], ratio[i]);
	}

	stop_timer(&timer);
	printf("Elapsed: %f seconds", get_timer_interval(&timer));

#ifndef LINUX
	system("pause");
#endif

	return 0;
}
