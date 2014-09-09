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
	double mse;
	struct cond_quantizer_list_t *quantizers;

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
	quantizers = generate_codebooks(&training_file, training_stats, dist, 0.5, BIT_ALLOC_MODE_NO_MIX, &mse);
	write_codebook(output_path, quantizers);

	stop_timer(&timer);
	printf("Expected MSE: %f\n", mse); 
	printf("Elapsed: %f seconds\n", get_timer_interval(&timer));

#ifndef LINUX
	system("pause");
#endif

	return 0;
}
