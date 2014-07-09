#ifndef _LINES_H_
#define _LINES_H_

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>

#include "pmf.h"

// This limits us to chunks that aren't too big to fit into a modest amount of memory at a time
#define MAX_LINES_PER_BLOCK			1000000
#define MAX_READS_PER_LINE			1022

// Error codes for reading a line block
#define LF_ERROR_NONE				0
#define LF_ERROR_NOT_FOUND			1
#define LF_ERROR_NO_MEMORY			2
#define LF_ERROR_TOO_LONG			4

/**
 * Points to a single line, which may be a pointer to a file in memory
 */
struct line_t {
	symbol_t *data;
};

/**
 * Points to a block of lines for incremental processing
 */
struct line_block_t {
	uint32_t count;
	struct line_t *lines;
};

/**
 * Points to a file descriptor that includes important metadata about the file
 */
struct quality_file_t {
	char *path;
	uint64_t lines;
	uint32_t columns;
	uint32_t block_count;
	struct line_block_t *blocks;
};

// Memory management
uint32_t load_file(const char *path, struct quality_file_t *info);

#endif
