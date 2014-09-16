/**
 * Utility functions for manipulating the data from files, like reading it into memory
 * and converting between the different formats we use
 */

#include "util.h"

#include <stdio.h>
#include <string.h>

#include "lines.h"

/**
 * This reads data from the given file pointer into memory, breaking it into segments
 * of the given number of lines, to ease memory management issues at the cost of some
 * overhead. This assumes that the file consists entirely of quality scores with no
 * other lines in between
 * @param path Path of the file to read
 * @param info Information structure to store in, this must be a valid pointer already
 * @param max_lines Maximum number of lines to read, will override the actual number in the file if >0
 */
uint32_t load_file(const char *path, struct quality_file_t *info, uint64_t max_lines) {
	uint32_t status, block_idx, line_idx;
	uint32_t j;
	char line[1024];
	FILE *fp;
	struct _stat finfo;

	// Load metadata into the info structure
	info->path = strdup(path);
	fp = fopen(path, "rt");
	if (!fp) {
		return LF_ERROR_NOT_FOUND;
	}

	// Use the first line to figure out how long the file is
	fgets(line, 1024, fp);
	info->columns = strlen(line) - 1;
	if (info->columns > MAX_READS_PER_LINE) {
		fclose(fp);
		return LF_ERROR_TOO_LONG;
	}

	// Figure out how many lines we'll need depending on whether we were limited or not
	_stat(path, &finfo);
	info->lines = finfo.st_size / ((uint64_t) (info->columns+1));
	if (max_lines > 0 && info->lines > max_lines) {
		info->lines = max_lines;
	}
	
	status = allocate_blocks(info);
	if (status != LF_ERROR_NONE)
		return status;

	// Process the file
	block_idx = 0;
	line_idx = 0;
	fseek(fp, 0, SEEK_SET);
	while (!feof(fp) && (block_idx * MAX_LINES_PER_BLOCK + line_idx) < info->lines) {
		// Read line and store in our array with data conversion, also stripping newlines
		fgets(line, 1024, fp);
		for (j = 0; j < info->columns; ++j) {
			info->blocks[block_idx].lines[line_idx].data[j] = line[j] - 33;
		}

		// Increment line/block pointers as necesary
		line_idx += 1;
		if (line_idx == info->blocks[block_idx].count) {
			line_idx = 0;
			block_idx += 1;
		}
	}

	fclose(fp);
	return LF_ERROR_NONE;
}

/**
 * Allocate an array of line block pointers and the memory within each block, so that we can
 * use it to store the results of reading the file
 */
uint32_t allocate_blocks(struct quality_file_t *info) {
	uint64_t lines_left = info->lines;
	symbol_t *sym_buf;
	struct line_block_t *cblock;
	struct line_t *cline;
	uint32_t i;

	// Figure out how many blocks we'll need to store this file
	info->block_count = (uint32_t) (info->lines / (uint64_t)MAX_LINES_PER_BLOCK);
	if (info->block_count * MAX_LINES_PER_BLOCK != info->lines) {
		info->block_count += 1;
	}

	info->blocks = (struct line_block_t *) calloc(info->block_count, sizeof(struct line_block_t));
	if (!info->blocks) {
		return LF_ERROR_NO_MEMORY;
	}
	cblock = info->blocks;

	while (lines_left > 0) {
		// Figure out how many lines we'll have in this block
		if (lines_left > MAX_LINES_PER_BLOCK) {
			lines_left -= MAX_LINES_PER_BLOCK;
			cblock->count = MAX_LINES_PER_BLOCK;
		}
		else {
			cblock->count = (uint32_t) lines_left;
			lines_left = 0;
		}

		// Allocate symbol buffer and array of line info structs for the block
		sym_buf = (symbol_t *) calloc(cblock->count, info->columns*sizeof(symbol_t));
		cblock->lines = (struct line_t *) calloc(cblock->count, sizeof(struct line_t));
		cline = cblock->lines;
		if (!sym_buf || !cblock->lines) {
			return LF_ERROR_NO_MEMORY;
		}

		// Save pointers into the symbol buffer for each line in the block
		for (i = 0; i < cblock->count; ++i) {
			cline->data = sym_buf;
			cline += 1;
			sym_buf += info->columns;
		}

		// Advance to the next line block
		cblock += 1;
	}

	return LF_ERROR_NONE;
}

/**
 * Deallocates the memory used to store file information in blocks
 */
void free_blocks(struct quality_file_t *info) {
	// Array of block pointers is a single allocation
	// For each block, array of lines and array of symbols are both single allocations

	uint32_t i;
	for (i = 0; i < info->block_count; ++i) {
		free(info->blocks[i].lines[0].data);
		free(info->blocks[i].lines);
	}
	free(info->blocks);
}
