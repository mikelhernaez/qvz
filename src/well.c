#include "well.h"

/**
 * Implementation of WELL-1024a PRNG
 * @param well_state_t RNG state to use for generation
 * @return 32 bit int random number
 */
uint32_t well_1024a(struct well_state_t *state) {
	uint32_t *s = state->state;
	uint32_t n = state->n;

	uint32_t z0 = s[(n+31)&31];
	uint32_t v_m1 = s[(n+3)&31];
	uint32_t v_m2 = s[(n+24)&31];
	uint32_t v_m3 = s[(n+10)&31];
	uint32_t z1 = s[n] ^ (v_m1 ^ (v_m1 >> 8));
	uint32_t z2 = (v_m2 ^ (v_m2 << 19)) ^ (v_m3 ^ (v_m3 << 14));

	s[n] = z1 ^ z2;
	n = (n + 31) & 31;
	s[n] = (z0 ^ (z0 << 11)) ^ (z1 ^ (z1 << 7)) ^ (z2 ^ (z2 << 13));
	state->n = n;
	return s[n];
}
