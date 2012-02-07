/*
 * This file is part of quip.
 *
 * Copyright (c) 2011 by Daniel C. Jones <dcjones@cs.washington.edu>
 *
 */


/*
 * twobit:
 * Nucleotide sequences encoded two bits per nucleotide.
 */

#ifndef QUIP_TWOBIT
#define QUIP_TWOBIT

#include "kmer.h"
#include <stdlib.h>
#include <stdio.h>

typedef struct twobit_t_ twobit_t;


twobit_t* twobit_alloc();
twobit_t* twobit_alloc_n(size_t n);
void      twobit_free(twobit_t*);
twobit_t* twobit_dup(const twobit_t*);
void      twobit_clear(twobit_t*);
void      twobit_reserve(twobit_t* s, size_t seqlen);

size_t twobit_len(const twobit_t*);
void   twobit_copy(twobit_t*, const char*);
void   twobit_copy_n(twobit_t*, const char*, size_t);
void   twobit_append(twobit_t*, const char*);
void   twobit_append_n(twobit_t*, const char*, size_t);
void   twobit_append_kmer(twobit_t*, kmer_t x, size_t k);
void   twobit_append_twobit(twobit_t*, const twobit_t*);
void   twobit_reverse(twobit_t*);


void   twobit_setc(twobit_t*, size_t i, char);
void   twobit_set(twobit_t*, size_t i, kmer_t);
kmer_t twobit_get(const twobit_t*, size_t i);
kmer_t twobit_get_kmer(const twobit_t*, size_t i, size_t k);
void   twobit_print(const twobit_t*, FILE*);
void   twobit_print_stdout(const twobit_t*);
int    twobit_cmp(const twobit_t*, const twobit_t*);
void   twobit_revcomp(twobit_t* dest, const twobit_t* src);

uint32_t twobit_hash(const twobit_t*);

#endif

