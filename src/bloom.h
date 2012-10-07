/*
 * This file is part of quip.
 *
 * Copyright (c) 2012 by Daniel C. Jones <dcjones@cs.washington.edu>
 *
 */

/*
 * bloom :
 * A variation of the "d-left counting bloom filter" proposed in:
 *
 *     Bonomi, F., Mitzenmacher, M., Panigrahy, R., Singh, S., & Varghese, G.
 *     (2006). An improved construction for counting bloom filters. 14th Annual
 *     European Symposium on Algorithms, LNCS 4168 (pp. 684–695). Springer.
 */


#ifndef QUIP_BLOOM
#define QUIP_BLOOM

#include "kmer.h"
#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>


typedef struct bloom_t_ bloom_t;

/* Allocate a new counting bloom filter, where n is the number of buckets per
 * table, and m is the number of cells per bucket.
 */
bloom_t* bloom_alloc(size_t n, size_t m);
void     bloom_clear(bloom_t*);
void     bloom_free(bloom_t*);

unsigned int bloom_inc(bloom_t*, kmer_t);
void         bloom_ldec(bloom_t*, kmer_t);
unsigned int bloom_add(bloom_t*, kmer_t, unsigned int d);
unsigned int bloom_get(bloom_t*, kmer_t);
void         bloom_del(bloom_t*, kmer_t);

#endif


/* The first thing we need to compress against a de bruijn graph is a bloom
 * filter that maps k-mers to conditional probabilities.
 *
 *
 * So, let's build a new data structure that maps k-mers to distributions.
 */



