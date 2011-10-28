
/* Represent a set over a fixed number of elements as a binary string */


#ifndef FAZER_BITSET_H
#define FAZER_BITSET_H

#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>

typedef struct bitset_t_ bitset_t;

bitset_t* bitset_alloc(size_t n);
void      bitset_free(bitset_t*);
void      bitset_copy(bitset_t*, const bitset_t*);
void      bitset_flip(bitset_t*, size_t i);
void      bitset_set(bitset_t*, size_t i, bool);
bool      bitset_get(const bitset_t*, size_t i);
void      bitset_on(bitset_t*, size_t i);
void      bitset_off(bitset_t*, size_t i);
void      bitset_allon(bitset_t*);
void      bitset_alloff(bitset_t*);
void      bitset_union(bitset_t*, const bitset_t*);
void      bitset_intersect(bitset_t*, const bitset_t*);
size_t    bitset_count_on(const bitset_t*);
void      bitset_print(const bitset_t*, FILE*);

#endif

