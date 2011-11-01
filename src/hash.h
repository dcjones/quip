/*
 * This file is part of fazer.
 *
 * Copyright (c) 2011 by Daniel C. Jones <dcjones@cs.washington.edu>
 *
 */

/*
 * hash :
 * A hash table for k-mers, borrowing ideas from:
 *
 *    Askitis, N., & Zobel, J. (2005). Cache-conscious collision resolution in
 *    string hash tables. String Processing and Information Retrieval (pp.
 *    91–102). Springer.
* 
 */


#ifndef FAZER_HASH
#define FAZER_HASH

#include "kmer.h"

typedef struct kmer_count_pair_t_
{
    kmer_t x;
    unsigned int count;
} kmer_count_pair_t;


typedef struct kmer_hash_t_ kmer_hash_t;

kmer_hash_t* kmer_hash_alloc();
void         kmer_hash_free(kmer_hash_t*);

size_t       kmer_hash_size(kmer_hash_t*);
void         kmer_hash_put(kmer_hash_t*, kmer_t, unsigned int);
unsigned int kmer_hash_get(kmer_hash_t*, kmer_t);

kmer_count_pair_t* kmer_hash_dump_sorted(kmer_hash_t*);

#endif


