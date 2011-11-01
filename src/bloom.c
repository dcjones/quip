
#include "bloom.h"
#include <assert.h>
#include <string.h>

/* TODO:
 * The big-endian code ought to be double checked, and tested.
 */


/* the number of subtables, hard-coded so I can use stack space in a few places
 * */
static const size_t d = 4;

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
    B->T = malloc(((d * n * m * cell_bytes) / 4 + 1) * 4);
    assert(B->T != NULL);
    memset(B->T, 0, ((d * n * m * cell_bytes) / 4 + 1) * 4);
    return B;
}


void bloom_free(bloom_t* B)
{
    free(B->T);
    free(B);
}



unsigned int bloom_get(bloom_t* B, kmer_t x)
{
    /* fingerprint */
    uint64_t h = kmer_hash(x);
#ifdef WORDS_BIGENDIAN
    uint32_t fp = (h & (uint64_t) fingerprint_mask);
#else
    uint32_t fp = (h & (uint64_t) fingerprint_mask);
#endif

    uint8_t* c;
    size_t i, j;
    for (i = 0; i < d; ++i) {
        h = kmer_hash_with_seed(x, h) % B->n;

        /* get bucket offset */
        c = B->T + (i * B->n + h) * B->m * cell_bytes;

        /* scan through cells */
        for (j = 0; j < B->m; ++j) {
            if (((*(uint32_t*) c) & fingerprint_mask) == fp) {
                return get_cell_count(c);
            }

            c += cell_bytes;
        }
    }

    return 0;
}

void bloom_inc(bloom_t* B, kmer_t x)
{
    /* fingerprint */
    uint64_t h = kmer_hash(x);
#ifdef WORDS_BIGENDIAN
    uint32_t fp = (h & (uint64_t) fingerprint_mask);
#else
    uint32_t fp = (h & (uint64_t) fingerprint_mask);
#endif
    uint32_t g;
    uint32_t cnt;

    uint8_t* cells[d];
    size_t bucket_sizes[d];

    uint8_t *c0, *c;
    size_t i, j;
    for (i = 0; i < d; ++i) {
        h = kmer_hash_with_seed(x, h) % B->n;

        /* get bucket offset */
        c0 = c = B->T + (i * B->n + h) * B->m * cell_bytes;

        /* scan through cells */
        for (j = 0; j < B->m; ++j) {
            g = (*(uint32_t*) c) & fingerprint_mask;

            assert(c < B->T + d * B->n * B->m * cell_bytes);

            if (g == fp) {
                cnt = get_cell_count(c);
                if (cnt < counter_mask) set_cell_count(c, cnt + 1);
                return;
            }
            else if (g == 0) {
                cells[i] = c;
                bucket_sizes[i] = (c - c0) / cell_bytes;
                break;
            }

            c += cell_bytes;
        }

        /* full bucket */
        if (j == B->m) {
            cells[i] = NULL;
            bucket_sizes[i] = B->m;
        }
    }

    /* find the smallest bucket, breaking ties to the left */
    size_t i_min = d;
    size_t min_bucket_size = B->m;
    for (i = 0; i < d && min_bucket_size > 0; ++i) {
        if (bucket_sizes[i] < min_bucket_size) {
            i_min = i;
            min_bucket_size = bucket_sizes[i];
        }
    }

    if (i_min < d) {
        (*(uint32_t*) cells[i_min]) = fp | 1; // figngerprint & count
    }
}


