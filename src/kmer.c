
#include "kmer.h"
#include "misc.h"
#include <stdlib.h>
#include <assert.h>


static const kmer_t complement1[4] = {3, 2, 1, 0};

/* reverse complements of 2, 4, and 8-kmers, resp. */
kmer_t* complement2 = NULL;
kmer_t* complement4 = NULL;
kmer_t* complement8 = NULL;


void kmer_init()
{
    kmer_t x;
    complement2 = malloc_or_die(0x100 * sizeof(kmer_t));
    for (x = 0; x <= 0xff; ++x) {
#if WORDS_BIGENDIAN
        complement2[x] = (complement1[x & 0x3] >> 2) | complement1[(x << 2)];
#else
        complement2[x] = (complement1[x & 0x3] << 2) | complement1[(x >> 2)];
#endif
    }


    complement4 = malloc_or_die(0x1000 * sizeof(kmer_t));
    for (x = 0; x <= 0xff; ++x) {
#if WORDS_BIGENDIAN
        complement4[x] = (complement2[x & 0xf] >> 4) | complement2[x << 4];
#else
        complement4[x] = (complement2[x & 0xf] << 4) | complement2[x >> 4];
#endif
    }


    complement8 = malloc_or_die(0x10000 * sizeof(kmer_t));
    for (x = 0; x < 0xffff; ++x) {
#if WORDS_BIGENDIAN
        complement8[x] = (complement4[x & 0xff] >> 8) | complement4[x << 8];
#else
        complement8[x] = (complement4[x & 0xff] << 8) | complement4[x >> 8];
#endif
    }
}


void kmer_free()
{
    free(complement2);
    free(complement4);
    free(complement8);
}


/* TODO:
 * mabye we can build a lookup table to do reverse complements
 * of very large chunks of the kmer at a time
 */


#if 0
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
#endif


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
        s[k - i - 1] = kmertochar((x & 0x3));
#ifdef WORDS_BIGENDIAN
        x <<= 2;
#else
        x >>= 2;
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
        y = (y >> 2) | complement1[x & 0x3];
        x <<= 2;
#else
        y = (y << 2) | complement1[x & 0x3]; 
        x >>= 2;
#endif
    }

    return y;
}


kmer_t kmer_revcomp(kmer_t x, size_t k)
{
    kmer_t y =
        (complement8[x & 0xffff] << 48) |
        (complement8[(x >> 16) & 0xffff] << 32) |
        (complement8[(x >> 32) & 0xffff] << 16) |
         complement8[x >> 48];

    return y >> (64 - (2 * k));
}


kmer_t kmer_canonical(kmer_t x, size_t k)
{
    kmer_t y = kmer_revcomp(x, k);
    return x < y ? x : y;
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
uint64_t kmer_hash_mix(uint64_t h1, uint64_t h2)
{
    const uint64_t c1 = 0x9ae16a3b2f90404fULL;
    const uint64_t c2 = 0x9ddfea08eb382d69ULL;

    h1 -= c1;
    uint64_t a = (h1 ^ h2) * c2;
    a ^= (a >> 47);
    uint64_t b = (h2 ^ a) * c2;
    b ^= (b >> 47);
    b *= c2;
    return b;
}


