
#include "bloom.h"
#include "misc.h"
#include <assert.h>
#include <string.h>

/* TODO:
 * The big-endian code ought to be double checked, and tested.
 */


/* the number of subtables, hard-coded so I can use stack space in a few places
 * */
#define NUM_SUBTABLES 4

/* careful, these numbers should not be changed independent of each other */
static const size_t   fingerprint_bits = 14;
static const uint32_t fingerprint_mask = 0xfffc00;
static const size_t   counter_bits     = 10;
static const uint32_t counter_mask     = 0x0003ff;
static const size_t   cell_bytes       = 3;


static uint32_t get_cell_count(uint8_t* c)
{
    return (*(uint32_t*) c) & counter_mask;
}


static void set_cell_count(uint8_t* c, uint32_t cnt)
{
    (*(uint32_t*) c) = ((*(uint32_t*) c) & fingerprint_mask) | (cnt & counter_mask);
}


struct bloom_t_
{
    uint8_t* T; 

    /* pointers into T, to save a little computation */
    uint8_t* subtable[NUM_SUBTABLES];

    /* number of buckets per subtable */
    size_t n;

    /* number of cells per bucket */
    size_t m;
};



bloom_t* bloom_alloc(size_t n, size_t m)
{
    bloom_t* B = malloc(sizeof(bloom_t));
    B->n = n;
    B->m = m;

    /* the '(. / 4 + 1) * 4' is to make sure things are aligned up to 32-bit
     * integers, mainly so valgrind doesn't whine. */
    B->T = malloc_or_die(((NUM_SUBTABLES * n * m * cell_bytes) / 4 + 1) * 4);
    memset(B->T, 0, ((NUM_SUBTABLES * n * m * cell_bytes) / 4 + 1) * 4);

    size_t i;
    for (i = 0; i < NUM_SUBTABLES; ++i) {
        B->subtable[i] = B->T + i * B->n * B->m * cell_bytes;
    }

    return B;
}


void bloom_clear(bloom_t* B)
{
    memset(B->T, 0, ((NUM_SUBTABLES * B->n * B->m * cell_bytes) / 4 + 1) * 4);
}


void bloom_free(bloom_t* B)
{
    free(B->T);
    free(B);
}



unsigned int bloom_get(bloom_t* B, kmer_t x)
{
    const size_t bytes_per_bucket = B->m * cell_bytes;

    /* fingerprint */
    uint64_t h, h1, h0 = kmer_hash(x);
    uint32_t fp = h0 & (uint64_t) fingerprint_mask;

    uint8_t* c;
    uint8_t* c_end;
    size_t i;
    h1 = h0;
    for (i = 0; i < NUM_SUBTABLES; ++i) {
        h1 = kmer_hash_mix(h0, h1);
        h = h1 % B->n;

        /* get bucket offset */
        c = B->subtable[i] + h * bytes_per_bucket;
        c_end = c + bytes_per_bucket;

        /* scan through cells */
        while (c < c_end) {
            if (((*(uint32_t*) c) & fingerprint_mask) == fp) {
                return get_cell_count(c);
            }

            c += cell_bytes;
        }
    }

    return 0;
}


void bloom_del(bloom_t* B, kmer_t x)
{
    const size_t bytes_per_bucket = B->m * cell_bytes;

    /* fingerprint */
    uint64_t h, h1, h0 = kmer_hash(x);
    uint32_t fp = h0 & (uint64_t) fingerprint_mask;

    uint8_t* c;
    uint8_t* c_end;
    size_t i;
    h1 = h0;
    for (i = 0; i < NUM_SUBTABLES; ++i) {
        h1 = kmer_hash_mix(h0, h1);
        h = h1 % B->n;

        /* get bucket offset */
        c = B->subtable[i] + h * bytes_per_bucket;
        c_end = c + bytes_per_bucket;

        /* scan through cells */
        while (c < c_end) {
            if (((*(uint32_t*) c) & fingerprint_mask) == fp) {
                (*(uint32_t*) c) &= ~(fingerprint_mask | counter_mask);
                return;
            }

            c += cell_bytes;
        }
    }
}


unsigned int bloom_inc(bloom_t* B, kmer_t x)
{
    return bloom_add(B, x, 1);
}


unsigned int bloom_add(bloom_t* B, kmer_t x, unsigned int d)
{
    const size_t bytes_per_bucket = B->m * cell_bytes;

    /* fingerprint */
    uint64_t h, h1, h0 = kmer_hash(x);
    uint32_t fp = h0 & (uint64_t) fingerprint_mask;

    uint32_t g;
    uint32_t cnt;

    uint8_t* cells[NUM_SUBTABLES];
    size_t bucket_sizes[NUM_SUBTABLES];

    uint8_t *c0, *c, *c_end;
    size_t i;
    h1 = h0;
    for (i = 0; i < NUM_SUBTABLES; ++i) {
        h1 = kmer_hash_mix(h0, h1);
        h = h1 % B->n;

        /* get bucket offset */
        c0 = c = B->subtable[i] + h * bytes_per_bucket;
        c_end = c + bytes_per_bucket;

        /* scan through cells */
        while (c < c_end) {
            g = (*(uint32_t*) c) & fingerprint_mask;

            if (g == fp) {
                cnt = get_cell_count(c);
                if (cnt + d < counter_mask) set_cell_count(c, cnt + d);
                else set_cell_count(c, counter_mask);
                return cnt + d;
            }
            else if (g == 0) {
                cells[i] = c;
                bucket_sizes[i] = (c - c0) / cell_bytes;
                break;
            }

            c += cell_bytes;
        }

        /* full bucket */
        if (c == c_end) {
            cells[i] = NULL;
            bucket_sizes[i] = B->m;
        }
    }

    /* find the smallest bucket, breaking ties to the left */
    size_t i_min = NUM_SUBTABLES;
    size_t min_bucket_size = B->m;
    for (i = 0; i < NUM_SUBTABLES && min_bucket_size > 0; ++i) {
        if (bucket_sizes[i] < min_bucket_size) {
            i_min = i;
            min_bucket_size = bucket_sizes[i];
        }
    }

    if (i_min < NUM_SUBTABLES) {
        if (d > counter_mask) d = counter_mask;
        (*(uint32_t*) cells[i_min]) = fp | d; // figngerprint & count
        return 1;
    }

    return 0; // no space for insertion
}


