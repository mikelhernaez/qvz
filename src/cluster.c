/**
 * k-means clustering implementation in C
 * 
 * The approach here is parallelizable and should be migrated to a block-based implementation
 * using opencl in order to run faster, or possibly just multithreaded, but I have left it
 * in plain C for the time being to get it working quickly. Note that the means established
 * are discrete values, rather than continuous.
 */

#include "util.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>

#define SYMBOLS 41
#define MAX_KMEANS_ITERATIONS 1000
#define MAX_COLUMN_SIZE 1024
#define MAX_LINES_PER_BLOCK 1000000

// Global constants, these should be moved into a metadata structure at some point
static uint32_t columns = 36;
static uint32_t total_lines = 0;

// Store a single line and the cluster to which it belongs
struct line_t {
	uint8_t cluster;			// Selected cluster for membership
	double *distances;			// Distance to each cluster center
	uint8_t *data;				// Data array for this line
};

// Stores information about each cluster
struct cluster_t {
	uint8_t id;					// Cluster ID
	uint8_t pad1;				// Alignment variable
	uint16_t pad2;				// Alignment variable

	uint32_t count;				// Number of lines in this cluster
	struct line_t **members;	// Array of pointers to the members

	struct line_t mean;			// Fake line whose values are the mean for this cluster
};

// Stores information about all clusters
struct cluster_list_t {
	uint8_t count;
	struct cluster_t *clusters;
};

// Stores information about a block of lines
struct line_block_t {
	uint32_t count;
	struct line_t *lines;
};

/**
 * Take a line and cluster information and calculates the distance, storing it in the line information vector
 */
void find_distance(struct line_t *line, struct cluster_t *cluster) {
	double d = 0.0;
	uint32_t i;
	uint32_t data, mean;

	for (i = 0; i < columns; ++i) {
		data = line->data[i];
		mean = cluster->mean.data[i];
		d += (data - mean) * (data - mean);
	}
	line->distances[cluster->id] = d;
}

/**
 * Assigns a cluster based on the one with the lowest distance
 */
uint8_t assign_cluster(struct line_t *line, struct cluster_list_t *clusters) {
	uint8_t id = line->cluster;
	uint8_t prev_id = line->cluster;
	uint8_t i;
	struct cluster_t *cluster;
	double d = line->distances[0];

	// Find the cluster with minimum distance
	for (i = 1; i < clusters->count; ++i) {
		if (line->distances[i] < d) {
			id = i;
			d = line->distances[i];
		}
	}

	// Assign to that cluster
	line->cluster = id;
	cluster = &clusters->clusters[id];
	cluster->members[cluster->count] = line;
	cluster->count += 1;

	return (prev_id == id) ? 0 : 1;
}

/**
 * Compare each line to each cluster to find distances
 */
uint8_t do_cluster_assignment(struct line_t *line, struct cluster_list_t *clusters) {
	uint8_t i;

	for (i = 0; i < clusters->count; ++i) {
		find_distance(line, &clusters->clusters[i]);
	}

	return assign_cluster(line, clusters);
}

/**
 * Calculates the mean for a cluster.
 */
void calculate_cluster_mean(struct cluster_t *cluster) {
	struct line_t *line;
	uint32_t i, j;
	uint64_t accumulator[MAX_COLUMN_SIZE];

	// Clear previous
	memset(accumulator, 0, columns*sizeof(uint64_t));

	// Sum up everything column-wise
	for (i = 0; i < cluster->count; ++i) {
		line = cluster->members[i];
		for (j = 0; j < columns; ++j) {
			accumulator[j] += line->data[j];
		}
	}

	// Integer division to find the mean, guaranteed to be less than SYMBOLS
	for (j = 0; j < columns; ++j) {
		cluster->mean.data[j] = (uint8_t) (accumulator[j] / cluster->count);
	}

}

/**
 * Updates the cluster means based on their assigned lines. Also clears the line count for
 * the next iteration.
 */
void recalculate_means(struct cluster_list_t *clusters) {
	uint8_t i;

	for (i = 0; i < clusters->count; ++i) {
		calculate_cluster_mean(&clusters->clusters[i]);
	}
}

/**
 * Calculates cluster assignments for the block of lines given, and return
 * status indicating that at least one line changed clusters
 */
uint8_t cluster_lines(struct line_block_t *block, struct cluster_list_t *clusters) {
	uint32_t i;
	uint8_t changed = 0;

	for (i = 0; i < block->count; ++i) {
		changed |= do_cluster_assignment(&block->lines[i], clusters);
	}

	return changed;
}

/**
 * Initialize the cluster means based on the data given, using random selection
 */
void initialize_kmeans_clustering(struct line_block_t *blocks, uint32_t block_count, struct cluster_list_t *clusters) {
	uint8_t j;
	uint32_t block_id;
	uint32_t line_id;

	for (j = 0; j < clusters->count; ++j) {
		block_id = rand() % block_count;
		line_id = rand() % blocks[block_id].count;
		memcpy(clusters->clusters[j].mean.data, blocks[block_id].lines[line_id].data, columns*sizeof(uint8_t));
		printf("Chose block %d, line %d.\n", block_id, line_id);
	}
}

/**
 * Allocate the memory used for the clusters based on the number wanted and column config
 */
void allocate_clusters(struct cluster_list_t *clusters, uint8_t count) {
	uint8_t j;

	// Allocate array of cluster structures
	clusters->count = count;
	clusters->clusters = (struct cluster_t *) calloc(count, sizeof(struct cluster_t));

	// Fill in each cluster
	for (j = 0; j < count; ++j) {
		clusters->clusters[j].id = j;
		clusters->clusters[j].count = 0;
		clusters->clusters[j].mean.data = (uint8_t *) calloc(columns, sizeof(uint8_t));
		// mean.cluster is skipped because it is unused
		// mean.distances remains a null pointer because it is unused
		clusters->clusters[j].members = (struct line_t **) calloc(total_lines, sizeof(struct line_t *));
	}
}

/**
 * Deallocate the memory used for the clusters. Only needed if we're not going to just exit the program
 * so it is a TODO for now.
 */
void deallocate_clusters(struct cluster_list_t *clusters) {}

/**
 * Allocate an array of line blocks able to store line pointers to the appropriate number of lines
 * along with the right kinds of distance buffer data
 */
struct line_block_t *allocate_line_block(uint32_t line_count, struct cluster_list_t *clusters, uint32_t *block_count) {
	uint32_t blocks;
	struct line_block_t *line_block;
	uint8_t *data_block;
	double *distance_block;
	uint32_t j, i;
	uint32_t lines;

	// Need ceiling block count
	blocks = (line_count / MAX_LINES_PER_BLOCK);
	if (blocks * MAX_LINES_PER_BLOCK < line_count)
		blocks += 1;

	// Allocate block array
	line_block = (struct line_block_t *) calloc(blocks, sizeof(struct line_block_t));
	*block_count = blocks;

	// For each block, allocate line structures, data and distance blocks, and configure the pointers
	for (j = 0; j < blocks; ++j) {
		if (line_count > MAX_LINES_PER_BLOCK)
			lines = MAX_LINES_PER_BLOCK;
		else
			lines = line_count;

		line_block[j].count = lines;
		line_block[j].lines = (struct line_t *) calloc(lines, sizeof(struct line_t));
		data_block = (uint8_t *) calloc(lines, columns * sizeof(uint8_t));
		distance_block = (double *) calloc(lines, clusters->count * sizeof(double));

		for (i = 0; i < lines; ++i) {
			line_block[j].lines[i].data = data_block;
			line_block[j].lines[i].distances = distance_block;
			data_block += columns;
			distance_block += clusters->count;
		}

		line_count -= lines;
	}

	return line_block;
}

/**
 * Deallocates an array of line blocks
 * TODO: Don't need to clean up yet!
 */
void deallocate_line_block(struct line_block_t *block) {}

/**
 * Do k-means clustering over the set of blocks given to produce a set of clusters that
 * fills the cluster list given
 */
void do_kmeans_clustering(struct line_block_t *blocks, uint32_t block_count, struct cluster_list_t *clusters) {
	uint32_t iter_count = 0;
	uint8_t changed = 1;
	uint32_t j;

	initialize_kmeans_clustering(blocks, block_count, clusters);

	while (iter_count < MAX_KMEANS_ITERATIONS && changed) {
		for (j = 0; j < clusters->count; ++j) {
			clusters->clusters[j].count = 0;
		}

		changed = 0;
		for (j = 0; j < block_count; ++j) {
			changed |= cluster_lines(&blocks[j], clusters);
		}

		recalculate_means(clusters);
		iter_count += 1;
		printf(".");
	}

	printf("\nTotal number of iterations: %d.\n", iter_count);
}

/**
 * Driver program that reads in the files and writes out the files
 */
int main(int argc, char **argv) {
	struct cluster_list_t clusters;
	struct line_block_t *blocks;
	struct line_t *qline;
	uint32_t block_count;
	FILE *fp;
	char line[1024];
	char filename[1024];
	struct _stat finfo;
	uint8_t cluster_count;
	uint32_t j, i, k;
	struct hrtimer_t timer;
	
	if (argc != 4) {
		printf("Usage: %s [file to cluster]  [number of clusters] [output prefix].\n", argv[0]);
		exit(1);
	}

	start_timer(&timer);
	srand(time(0));

	fp = fopen(argv[1], "rt");
	if (!fp) {
		perror("Unable to open input file");
		exit(1);
	}

	// Figure out how many columns we have from the first line
	fgets(line, 1024, fp);
	columns = strlen(line)-1;

	// Then find the number of lines from the file size
	_stat(argv[1], &finfo);
	total_lines = finfo.st_size / (columns+1);

	// Initialize cluster and line block storage
	cluster_count = (uint8_t) atoi(argv[2]);
	allocate_clusters(&clusters, cluster_count);
	blocks = allocate_line_block(total_lines, &clusters, &block_count);

	// Read in the file to our line blocks
	j = 0;
	i = 0;
	do {
		memcpy(blocks[j].lines[i].data, line, columns*sizeof(uint8_t));

		// Adjust quality scores by 33
		for (k = 0; k < columns; ++k) {
			blocks[j].lines[i].data[k] -= 33;
		}

		// Find the place where we're storing the next line
		i += 1;
		if (i == blocks[j].count) {
			j += 1;
			i = 0;
		}

		// Read in this line
		fgets(line, 1024, fp);
	} while (!feof(fp));
	fclose(fp);

	// Run k-means clustering
	printf("Beginning clustering...\n");
	do_kmeans_clustering(blocks, block_count, &clusters);
	printf("Clustering done.\n");

	// Write out each cluster to a file
	for (j = 0; j < clusters.count; ++j) {
		sprintf(filename, "%s_cluster_%d.txt", argv[3], j);
		printf("Cluster %d: ", j);
		for (k = 0; k < columns; ++k) {
			clusters.clusters[j].mean.data[k] += 33;
			printf("%c", clusters.clusters[j].mean.data[k]);
		}
		printf("\n");
		printf("Members: %d.\n", clusters.clusters[j].count);
		
		fp = fopen(filename, "wt");
		for (i = 0; i < clusters.clusters[j].count; ++i) {
			qline = clusters.clusters[j].members[i];

			// Adjust quality scores for this line
			for (k = 0; k < columns; ++k) {
				qline->data[k] += 33;
			}

			// Write the line with a newline
			fwrite(qline->data, columns, sizeof(uint8_t), fp);
			fwrite("\n", 1, sizeof(char), fp);
		}
		fclose(fp);
		printf("Write cluster %d to file %s.\n", j, filename);
	}

	stop_timer(&timer);
	printf("Time elapsed: %f sec.\n", get_timer_interval(&timer));

	return 0;
}
