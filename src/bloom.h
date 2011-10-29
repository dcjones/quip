
/*
 * A variation of the "d-left counting bloom filter" proposed by Bonomi, et al.
 * in the paper entitled: "An Improved Construction for Counting Bloom Filters".
 *
 */


#ifndef FAZER_BLOOM
#define FAZER_BLOOM

#include "kmer.h"
#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>


typedef struct bloom_t_ bloom_t;
typedef uint8_t bloom_count_t;

/* Allocate a new counting bloom filter, where n is the number of buckets per
 * table, and m is the number of cells per bucket.
 */
bloom_t* bloom_alloc(size_t n, size_t m);
void     bloom_free(bloom_t*);

void         bloom_inc(bloom_t*, kmer_t);
unsigned int bloom_get(bloom_t*, kmer_t);

#endif


