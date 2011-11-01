
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


static kmer_t kmer_high_order(kmer_t x, size_t k)
{
#ifdef WORDS_BIGENDIAN
    return x << (2 * (k - 1));
#else
    return x >> (2 * (k - 1));
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

char kmertochar(kmer_t x)
{
    switch (x & 0x3) {
        case 0: return 'A';
        case 1: return 'C';
        case 2: return 'G';
        case 3: default: return 'T';
    }
}


kmer_t strtokmer(const char* s)
{
    size_t k = 0;
    kmer_t x_i, x = 0;
    while (*s && k++ < 4 * sizeof(kmer_t)) {
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


void kmertostr(kmer_t x, char* s, size_t k)
{
    size_t i = 0;
    for (i = 0; i < k; ++i) {
        s[i] = kmertochar((x & 0x3));
#ifdef WORDS_BIGENDIAN
        x >>= 2;
#else
        x <<= 2;
#endif
    }

    s[i] = '\0';
}


kmer_t kmer_get_nt(kmer_t* x, size_t i)
{
    size_t idx = (4 * sizeof(kmer_t)) / i;
    size_t off = (4 * sizeof(kmer_t)) % i;

#ifdef WORDS_BIGENDIAN
    return 0x3 & (x[idx] << off);
#else
    return 0x3 & (x[idx] >> off);
#endif
}


kmer_t kmer_comp(kmer_t x, size_t k)
{
    kmer_t y = 0;
    while (k--) {
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


kmer_t kmer_revcomp(kmer_t x, size_t k)
{
    kmer_t y = 0;
    while (k--) {
#ifdef WORDS_BIGENDIAN
        y = (y >> 2) | kmer_comp1(x << (2 * (k - 1)));
        x >>= 2;
#else
        y = (y << 2) | kmer_comp1(x >> (2 * (k - 1)));
        x <<= 2;
#endif
    }

    return y;
}


kmer_t kmer_canonical(kmer_t x, size_t k)
{
    /* to decide whether x or revcomp(x) is lexigraphically smaller, we only
     * need to look at the higher order nucleotide and complement of the low
     * order nucleotide. */

    if (kmer_high_order(x, k) < kmer_comp1(kmer_low_order(x))) {
        return x;
    }
    else {
        return kmer_revcomp(x, k);
    }
}



/* This is Thomas Wang's hash function for 64-bit integers. */
uint64_t kmer_hash(kmer_t x)
{
    x = (~x) + (x << 21);
    x = x ^ (x >> 24);
    x = (x + (x << 3)) + (x << 8); // x * 265
    x = x ^ (x >> 14);
    x = (x + (x << 2)) + (x << 4); // x * 21
    x = x ^ (x >> 28);
    x = x + (x << 31);
    return x;
}


/* This is taken from the Hash128to64 function in CityHash */
uint64_t kmer_hash_with_seed(kmer_t x, uint64_t h2)
{
    const uint64_t c1 = 0x9ae16a3b2f90404fULL;
    const uint64_t c2 = 0x9ddfea08eb382d69ULL;

    uint64_t h1 = kmer_hash(x) - c1;
    uint64_t a = (h1 ^ h2) * c2;
    a ^= (a >> 47);
    uint64_t b = (h2 ^ a) * c2;
    b ^= (b >> 47);
    b *= c2;
    return b;
}


