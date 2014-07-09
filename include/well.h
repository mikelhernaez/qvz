#ifndef _WELL_H_
#define _WELL_H_
/**
 * Declarations necessary for using the WELL-1024a PRNG in other code
 */

#include <stdint.h>

struct well_state_t {
	uint32_t state[32];
	uint32_t n;
};


uint32_t well_1024a(struct well_state_t *state);

#endif
