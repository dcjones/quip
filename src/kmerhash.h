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



typedef struct kmerhash_t_ kmerhash_t;

typedef struct kmer_pos_t_
{
    uint32_t contig_idx;
    int32_t contig_pos;
} kmer_pos_t;

kmerhash_t* kmerhash_alloc();
void         kmerhash_free(kmerhash_t*);

size_t kmerhash_size(kmerhash_t*);
void   kmerhash_put(kmerhash_t*, kmer_t, uint32_t contig_idx, int32_t contig_pos);
size_t kmerhash_get(kmerhash_t*, kmer_t, kmer_pos_t**);

#endif


