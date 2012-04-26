/*
 * This file is part of quip.
 *
 * Copyright (c) 2012 by Daniel C. Jones <dcjones@cs.washington.edu>
 *
 */


/*
 * seqset:
 * Fixed size set of nucleotide sequences with associated counts.
 */


#ifndef QUIP_SEQSET
#define QUIP_SEQSET

#include "twobit.h"

typedef struct seqset_t_ seqset_t;

typedef struct seqset_value_t_
{
    twobit_t* seq;
    uint32_t  cnt;
} seqset_value_t;


seqset_t* seqset_alloc();
seqset_t* seqset_alloc_fixed(size_t N);
void      seqset_clear(seqset_t*);
void      seqset_free(seqset_t*);
uint32_t  seqset_inc(seqset_t*, const twobit_t*);
size_t    seqset_size(const seqset_t*);

seqset_value_t* seqset_dump(const seqset_t*);

typedef struct
{
    const seqset_value_t* val;

    const seqset_t* S;
    size_t i;
} seqset_iter_t;


seqset_iter_t* seqset_iter_alloc();
void           seqset_iter_free(seqset_iter_t*);

void seqset_iter_start(seqset_iter_t*, const seqset_t*);
void seqset_iter_next(seqset_iter_t*);
bool seqset_iter_finished(const seqset_iter_t*);


#endif

