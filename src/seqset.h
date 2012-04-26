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

typedef struct
{
    /* We store every sequence we can in two bit per nucelotide, but those with
     * N's are kept as eight bits per nucleotide. */
    union {
        twobit_t* tb;
        char*     eb;
    } seq;

    bool is_twobit;

    /* Number of accurances of the given sequence. */
    uint32_t  cnt;

    /* We assign indexes sequentially as new sequences are inserted, but they
     * are not robust to deletion. */
    uint32_t  idx;
} seqset_value_t;


seqset_t* seqset_alloc();
void      seqset_clear(seqset_t*);
void      seqset_free(seqset_t*);
uint32_t  seqset_inc_tb(seqset_t*, const twobit_t*);
uint32_t  seqset_inc_eb(seqset_t*, const char*);
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

