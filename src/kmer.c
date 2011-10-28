
#include "kmer.h"
#include <stdlib.h>


static kmer_t kmer_comp1(kmer_t x)
{
    switch(x & 0x3) {
        case 0: return 3;
        case 1: return 2;
        case 2: return 1;
        case 3: default: return 0;
    }
}


static kmer_t kmer_high_order(kmer_t x)
{
#ifdef WORDS_BIGENDIAN
    return x << (sizeof(kmer_t) - 2);
#else
    return x >> (sizeof(kmer_t) - 2);
#endif
}


static kmer_t kmer_low_order(kmer_t x)
{
    return x & 0x3;
}


kmer_t chartokmer(char c)
{
    switch (c) {
        case 'a': case 'A': case '0': return 0;
        case 'c': case 'C': case '1': return 1;
        case 'g': case 'G': case '2': return 2;
        case 't': case 'T': case '3': return 3;

        /* error reporting is a bit subtle, remember to check this */
        default: return 4;
    }
}


kmer_t strtokmer(const char* s)
{
    kmer_t x_i, x = 0;
    while (*s) {
        x_i = chartokmer(*s);
        if (x_i > 3) break;
#ifdef WORDS_BIGENDIAN
        x = (x >> 2) | x_i;
#else
        x = (x << 2) | x_i;
#endif
    }

    return x;
}


kmer_t kmer_comp(kmer_t x)
{
    kmer_t y = 0;
    size_t n = 4 * sizeof(kmer_t);
    while (n--) {
#ifdef WORDS_BIGENDIAN
        y = (y >> 2) | kmer_comp1(x);
        x <<= 2;
#else
        y = (y << 2) | kmer_comp1(x); 
        x >>= 2;
#endif
    }

    return y;
}


kmer_t kmer_revcomp(kmer_t x)
{
    kmer_t y = 0;
    size_t n = 4 * sizeof(kmer_t);
    while (n--) {
#ifdef WORDS_BIGENDIAN
        y = (y >> 2) | kmer_comp1((x << (sizeof(kmer_t) - 2)) & 0x3);
        x >>= 2;
#else
        y = (y << 2) | kmer_comp1((x >> (sizeof(kmer_t) - 2)) & 0x3);
        x <<= 2;
#endif
    }

    return y;
}


kmer_t kmer_canonical(kmer_t x)
{
    /* to decide whether x or revcomp(x) is lexigraphically smaller, we only
     * need to look at the higher order nucleotide and complement of the low
     * order nucleotide. */

    if (kmer_high_order(x) < kmer_comp1(kmer_low_order(x))) {
        return x;
    }
    else {
        return kmer_revcomp(x);
    }
}




