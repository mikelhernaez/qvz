/**
 * k-means clustering implementation in C
 * 
 * The approach here is parallelizable and should be migrated to a block-based implementation
 * using opencl in order to run faster, or possibly just multithreaded, but I have left it
 * in plain C for the time being to get it working quickly. Note that the means established
 * are discrete values, rather than continuous.
 */

#include "util.h"

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <malloc.h>

#include "pmf.h"
#include "codebook.h"
#include "cluster.h"

/**
 * Allocate the memory used for the clusters based on the number wanted and column config
 */
struct cluster_list_t *alloc_cluster_list(struct quality_file_t *info) {
	uint8_t j;
	struct cluster_list_t *rtn = (struct cluster_list_t *) calloc(1, sizeof(struct cluster_list_t));

	// Allocate array of cluster structures
	rtn->count = info->cluster_count;
	rtn->clusters = (struct cluster_t *) calloc(info->cluster_count, sizeof(struct cluster_t));

	// Fill in each cluster
	for (j = 0; j < info->cluster_count; ++j) {
		rtn->clusters[j].id = j;
		rtn->clusters[j].count = 0;
		rtn->clusters[j].mean.data = (symbol_t *) calloc(info->columns, sizeof(symbol_t));
		// mean.cluster is skipped because it is unused
		// mean.distances remains a null pointer because it is unused
		rtn->clusters[j].members = (struct line_t **) calloc(info->lines, sizeof(struct line_t *));
		rtn->clusters[j].training_stats = alloc_conditional_pmf_list(info->alphabet, info->columns);
	}

	return rtn;
}

/**
 * Deallocate the memory used for the clusters.
 */
void free_cluster_list(struct cluster_list_t *clusters) {
	uint8_t j;

	for (j = 0; j < clusters->count; ++j) {
		free(clusters->clusters[j].mean.data);
		free(clusters->clusters[j].members);
		free_conditional_pmf_list(clusters->clusters[j].training_stats);
	}
	free(clusters->clusters);
	free(clusters);
}

/**
 * Calculates cluster assignments for the block of lines given, and return
 * status indicating that at least one line changed clusters
 */
uint8_t cluster_lines(struct line_block_t *block, struct quality_file_t *info) {
	uint32_t i;
	uint8_t changed = 0;

	for (i = 0; i < block->count; ++i) {
		changed |= do_cluster_assignment(&block->lines[i], info);
	}

	return changed;
}

/**
 * Updates the cluster means based on their assigned lines. Also clears the line count for
 * the next iteration.
 */
void recalculate_means(struct quality_file_t *info) {
	uint8_t i;

	for (i = 0; i < info->cluster_count; ++i) {
		calculate_cluster_mean(&info->clusters->clusters[i], info);
	}
}

/**
 * Calculates the mean for a cluster.
 */
void calculate_cluster_mean(struct cluster_t *cluster, struct quality_file_t *info) {
	struct line_t *line;
	uint32_t i, j;
	uint64_t *accumulator = _alloca(info->columns*sizeof(uint64_t));

	// Clear previous
	memset(accumulator, 0, info->columns*sizeof(uint64_t));

	// Sum up everything column-wise
	for (i = 0; i < cluster->count; ++i) {
		line = cluster->members[i];
		for (j = 0; j < info->columns; ++j) {
			accumulator[j] += line->data[j];
		}
	}

	// Integer division to find the mean, guaranteed to be less than the alphabet size
	for (j = 0; j < info->columns; ++j) {
		cluster->mean.data[j] = (uint8_t) (accumulator[j] / cluster->count);
	}
}

/**
 * Compare each line to each cluster to find distances
 */
uint8_t do_cluster_assignment(struct line_t *line, struct quality_file_t *info) {
	uint8_t i;

	for (i = 0; i < info->cluster_count; ++i) {
		find_distance(line, &info->clusters->clusters[i], info);
	}

	return assign_cluster(line, info);
}

/**
 * Assigns a cluster based on the one with the lowest distance
 */
uint8_t assign_cluster(struct line_t *line, struct quality_file_t *info) {
	uint8_t id = 0;
	uint8_t prev_id = line->cluster;
	uint8_t i;
	struct cluster_t *cluster;
	double d = line->distances[0];

	// Find the cluster with minimum distance
	for (i = 1; i < info->cluster_count; ++i) {
		if (line->distances[i] < d) {
			id = i;
			d = line->distances[i];
		}
	}

	// Assign to that cluster
	line->cluster = id;
	cluster = &info->clusters->clusters[id];
	cluster->members[cluster->count] = line;
	cluster->count += 1;

	return (prev_id == id) ? 0 : 1;
}

/**
 * Take a line and cluster information and calculates the distance, storing it in the line information vector
 */
void find_distance(struct line_t *line, struct cluster_t *cluster, struct quality_file_t *info) {
	double d = 0.0;
	uint32_t i;
	uint32_t data, mean;

	for (i = 0; i < info->columns; ++i) {
		data = line->data[i];
		mean = cluster->mean.data[i];
		//d += (data - mean) * (data - mean);
		d += get_distortion(info->dist, data, mean);
	}
	line->distances[cluster->id] = d;
}

/**
 * Initialize the cluster means based on the data given, using random selection
 */
void initialize_kmeans_clustering(struct quality_file_t *info) {
	uint8_t j;
	uint32_t block_id;
	uint32_t line_id;
	struct cluster_list_t *clusters = info->clusters;

	for (j = 0; j < info->cluster_count; ++j) {
		block_id = rand() % info->block_count;
		line_id = rand() % info->blocks[block_id].count;
		memcpy(clusters->clusters[j].mean.data, info->blocks[block_id].lines[line_id].data, info->columns*sizeof(uint8_t));
		if (info->opts->verbose) {
			printf("Chose block %d, line %d.\n", block_id, line_id);
		}
	}
}

/**
 * Do k-means clustering over the set of blocks given to produce a set of clusters that
 * fills the cluster list given
 */
void do_kmeans_clustering(struct quality_file_t *info) {
	uint32_t iter_count = 0;
	uint8_t changed = 1;
	uint32_t j;
	struct cluster_list_t *clusters = info->clusters;

	initialize_kmeans_clustering(info);

	while (iter_count < MAX_KMEANS_ITERATIONS && changed) {
		for (j = 0; j < clusters->count; ++j) {
			clusters->clusters[j].count = 0;
		}

		changed = 0;
		for (j = 0; j < info->block_count; ++j) {
			changed |= cluster_lines(&info->blocks[j], info);
		}

		recalculate_means(info);
		iter_count += 1;
		if (info->opts->verbose) {
			printf(".");
		}
	}

	if (info->opts->verbose) {
		printf("\nTotal number of iterations: %d.\n", iter_count);
	}
}
