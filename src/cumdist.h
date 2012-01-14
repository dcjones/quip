/*
 * This file is part of quip.
 *
 * Copyright (c) 2011 by Daniel C. Jones <dcjones@cs.washington.edu>
 *
 */

/*
 * cumdist:
 * A representation of a cumulative probability distribution that can be
 * efficiently updated.
 *
 * I think this is originally due to:
 *     P.M. Fenwick, "A new data structure for cumulative frequency tables,"
 *     Softw. Pract. Exper., vol. 24, pp. 327--336, March 1994
 *
 */

#ifndef QUIP_CUMDIST
#define QUIP_CUMDIST

#include <stdlib.h>
#include <stdint.h>

typedef struct cumdist_t_ cumdist_t;


cumdist_t* cumdist_alloc(size_t n);
void       cumdist_free(cumdist_t*);

/* frequency */
uint32_t cumdist_p(const cumdist_t*, size_t i);
uint32_t cumdist_p_norm(const cumdist_t*, size_t i);

/* the frequency of elements <= i, where the ordering is fixed but arbitrary
 * (i.e., DO NOT depend on i < j implying P(i) <= P(j))
 */
uint32_t cumdist_P(const cumdist_t*, size_t i);
uint32_t cumdist_P_norm(const cumdist_t*, size_t i);

/* total frequency */
uint32_t cumdist_Z(const cumdist_t*);

/* add x to thefrequency of element i */
void cumdist_add(cumdist_t*, size_t i, uint32_t x);


#endif

