/*
 * This file is part of quip.
 *
 * Copyright (c) 2012 by Daniel C. Jones <dcjones@cs.washington.edu>
 *
 */

/*
 * dist:
 * Efficient representation of discrete probability distributions for arithmetic
 * coding.
 */


#ifndef QUIP_DIST
#define QUIP_DIST

#include <stdlib.h>
#include <stdint.h>

extern const size_t dist_length_shift;

typedef struct dist_t_
{
    /* alphabet size */
    size_t n;

    /* cumulative symbol frequency */
    uint32_t* ps;

    /* symbol counts */
    uint32_t* cs;

    /* total number of occurances */
    uint32_t z;

    /* decoder table */
    uint32_t* dec;
    size_t dec_size;
    size_t dec_shift;

    /* number of new observations until the distribution is updated */
    size_t update_delay;
} dist_t;


/* Integer used to represent symbols in the alphabet. */
typedef uint32_t symb_t;

/* Allocate a distribution over the alphabet [0, n - 1].
 * The structure is either allocated with the specify intent of either decoding
 * or encoding to save memory.
 * */
dist_t* dist_alloc_encode(size_t n);
dist_t* dist_alloc_decode(size_t n);
void    dist_free(dist_t*);

/* alphabet size */
size_t dist_n(const dist_t*);

/* cumulative symbol frequency */
uint32_t dist_P(const dist_t*, symb_t x);

/* add k to the number of occurances of symbol x */
void dist_add(dist_t*, symb_t x, uint32_t k);

/* update distribution to reflect calls to dist_add */
void dist_update(dist_t* D);

/* find the symbol corresponding to the given frequency */
size_t dist_find(dist_t*, uint32_t);

#endif

