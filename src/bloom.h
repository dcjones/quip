
#ifndef FAZER_BLOOM
#define FAZER_BLOOM

#include "kmer.h"
#include <stdlib.h>

typedef struct bloom_t_ bloom_t;

bloom_t* bloom_alloc(size_t m, size_t k);
void     bloom_free(bloom_t*);
void     bloom_get(bloom_t*, kmer_t);

#endif


