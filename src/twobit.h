/*
 * This file is part of fazer.
 *
 * Copyright (c) 2011 by Daniel C. Jones <dcjones@cs.washington.edu>
 *
 */


/*
 * twobit:
 * Nucleotide sequences encoded two bits per nucleotide.
 */

#ifndef FAZER_TWOBIT
#define FAZER_TWOBIT

#include "kmer.h"
#include <stdlib.h>

typedef struct twobit_t_ twobit_t;


twobit_t* twobit_alloc();
twobit_t* twobit_alloc_n(size_t n);
void      twobit_free(twobit_t*);


size_t twobit_len(twobit_t*);
void   twobit_append(twobit_t*, const char*);
void   twobit_append_n(twobit_t*, const char*, size_t);


void   twobit_set(size_t i, twobit_t*, char);
kmer_t twobit_get(size_t i, twobit_t*);



#endif

