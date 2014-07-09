#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "pmf.h"
#include "distortion.h"
#include "quantizer.h"

int main(int argc, char **argv) {
	struct distortion_t *dist = generate_distortion_matrix(41, DISTORTION_MSE);
	struct alphabet_t *alphabet = alloc_alphabet(41);
	struct pmf_t *pmf = alloc_pmf(alphabet);
	struct quantizer_t *q;
	uint8_t i;
	double sum = 0.0;
	double t;

	for (i = 0; i < 41; ++i) {
		t = (i - 20)/3.;
		pmf->pmf[i] = exp(-t*t);
		sum += pmf->pmf[i];
	}
	for (i = 0; i < 41; ++i) {
		pmf->pmf[i] = pmf->pmf[i] / sum;
	}
	pmf->pmf_ready = 1;

	print_alphabet(alphabet);
	print_pmf(pmf);

	q = generate_quantizer(pmf, dist, 4, &t);

	print_quantizer(q);
	printf("Expected distortion: %f.\n", t);

#ifndef LINUX
	system("pause");
#endif

	return 0;
}
