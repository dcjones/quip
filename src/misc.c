/*
 * This file is part of quip.
 *
 * Copyright (c) 2012 by Daniel C. Jones <dcjones@cs.washington.edu>
 *
 */


#include "misc.h"
#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>



void* malloc_or_die(size_t n)
{
    void* p = malloc(n);
    if (p == NULL) {
        fprintf(stderr, "Can not allocate %zu bytes.\n", n);
        exit(EXIT_FAILURE);
    }
    return p;
}


void* realloc_or_die(void* ptr, size_t n)
{
    void* p = realloc(ptr, n);
    if (p == NULL) {
        fprintf(stderr, "Can not allocate %zu bytes.\n", n);
        exit(EXIT_FAILURE);
    }
    return p;
}


FILE* fopen_or_die(const char* path, const char* mode)
{
    FILE* f = fopen(path, mode);
    if (f == NULL) {
        fprintf(stderr, "Can not open file %s with mode %s.\n", path, mode);
        exit(EXIT_FAILURE);
    }
    return f;
}


/* This is MurmurHash3. The original C++ code was placed in the public domain
 * by its author, Austin Appleby. */

static inline uint32_t fmix(uint32_t h)
{
    h ^= h >> 16;
    h *= 0x85ebca6b;
    h ^= h >> 13;
    h *= 0xc2b2ae35;
    h ^= h >> 16;

    return h;
}


static inline uint32_t rotl32(uint32_t x, int8_t r)
{
    return (x << r) | (x >> (32 - r));
}


uint32_t murmurhash3(const uint8_t* data, size_t len_)
{
    const int len = (int) len_;
    const int nblocks = len / 4;

    uint32_t h1 = 0xc062fb4a;

    uint32_t c1 = 0xcc9e2d51;
    uint32_t c2 = 0x1b873593;

    //----------
    // body

    const uint32_t * blocks = (const uint32_t*) (data + nblocks * 4);

    int i;
    for(i = -nblocks; i; i++)
    {
        uint32_t k1 = blocks[i];

        k1 *= c1;
        k1 = rotl32(k1, 15);
        k1 *= c2;

        h1 ^= k1;
        h1 = rotl32(h1, 13); 
        h1 = h1*5+0xe6546b64;
    }

    //----------
    // tail

    const uint8_t * tail = (const uint8_t*)(data + nblocks*4);

    uint32_t k1 = 0;

    switch(len & 3)
    {
        case 3: k1 ^= tail[2] << 16;
        case 2: k1 ^= tail[1] << 8;
        case 1: k1 ^= tail[0];
              k1 *= c1; k1 = rotl32(k1,15); k1 *= c2; h1 ^= k1;
    }

    //----------
    // finalization

    h1 ^= len;

    h1 = fmix(h1);

    return h1;
}


static unsigned char complement(unsigned char c)
{
    switch (c) {
        case 'a': return 't';
        case 'c': return 'g';
        case 'g': return 'c';
        case 't': return 'a';
        case 'n': return 'n';
        case 'A': return 'T';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'T': return 'A';
        case 'N': return 'N';
        default:  return c;
    }
}


void str_revcomp(unsigned char* seq, size_t n)
{
    char c;
    size_t i, j;
    i = 0;
    j = n - 1;
    while (i < j) {
        c = complement(seq[i]);
        seq[i] = complement(seq[j]);
        seq[j] = c;
        i++; j--;
    }

    if (i == j) seq[i] = complement(seq[i]);
}


void str_rev(unsigned char* seq, size_t n)
{
    char c;
    size_t i, j;
    i = 0;
    j = n - 1;
    while (i < j) {
        c = seq[i];
        seq[i] = seq[j];
        seq[j] = c;
        i++; j--;
    }
}

#ifndef HAVE_VASPRINTF

int vasprintf(char **ret, const char *format, va_list args)
{
    va_list copy;
    va_copy(copy, args);

    /* Make sure it is determinate, despite manuals indicating otherwise */
    *ret = 0;

    int count = vsnprintf(NULL, 0, format, args);
    if (count >= 0)
    {
        char* buffer = malloc(count + 1);
        if (buffer != NULL)
        {
            count = vsnprintf(buffer, count + 1, format, copy);
            if (count < 0)
                free(buffer);
            else
                *ret = buffer;
        }
    }
    va_end(copy);  // Each va_start() or va_copy() needs a va_end()

    return count;
}

#endif

