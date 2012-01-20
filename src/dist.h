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

typedef struct dist_t_ dist_t;

/* allocate a distribution over [0, n - 1] */
dist_t* dist_alloc(size_t n);
void    dist_free(dist_t*);

/* symbol frequency */
uint32_t dist_P(const dist_t*, size_t k);

/* cumulative symbol frequency */
uint32_t dist_C(const dist_t*, size_t k);

/* add x to the number of occurances of symbol k */
void dist_add(dist_t*, size_t k, uint32_t x);

/* apply changes made from one or more calls to dist_add */
void dist_update(dist_t*);

/* find the symbol corresponding to the given frequency */
size_t dist_find(dist_t*, uint32_t);

#endif

