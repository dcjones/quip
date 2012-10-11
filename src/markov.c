
#include <math.h>
#include <string.h>

#include "dist.h"
#include "markov.h"
#include "misc.h"

/* Number of sub-tables used by the dlCFB. */
#define NUM_SUBTABLES 4

/* Number of cells per bucket. */
#define NUM_CELLS_PER_BUCKET 8


/* Each cell holds a 16-bit fingerprint (hash) of the k-mer as well as a
 * conditional distribution over nucleotides following that k-mer. */
typedef struct cell_t_
{
    uint16_t fingerprint;
    dist16_t dist;
} cell_t;


/* Compute the fingerprint from a hash of a k-mer.
 *
 * We compute a regular k-mer hash, take the low order 16-bits and replace 0
 * with 1, since we use 0 to represent an unused cell. */
static uint16_t fingerprint_hash(uint64_t h)
{
    h &= 0xffff;
    return h == 0 ? 1 : h;
}


/* The idea here is to use the dlCBF to map a k-mer to a distribution. */
struct markov_t_
{
    cell_t* subtables[NUM_SUBTABLES];

    /* If a k-mer could not be found or inserted in the dlcbf, encode against
     * a lower-order dense markov chain. */
    cond_dist16_t catchall;

    /* Order of the catchall markov chain. */
    size_t k_catchall;

    /* Number of buckets per sub-table. */
    size_t n;

    /* Order of the markov chain */
    size_t k;

    /* Bitmask for the (k-1)-mer  */
    kmer_t xmask;
    kmer_t xmask_catchall;
};


markov_t* markov_create(size_t n, size_t k, size_t k_catchall)
{
    markov_t* mc = malloc_or_die(sizeof(markov_t));
    mc->n = (size_t) ceil((double) n / (double) NUM_CELLS_PER_BUCKET
                                     / (double) NUM_SUBTABLES);
    mc->k = k;
    mc->xmask = kmer_mask(k);
    mc->k_catchall = k_catchall;
    mc->xmask_catchall = kmer_mask(k_catchall);

    cond_dist16_init(&mc->catchall, 1 << (2 * k_catchall));

    size_t i;
    size_t N = mc->n * NUM_CELLS_PER_BUCKET * sizeof(cell_t);
    for (i = 0; i < NUM_SUBTABLES; ++i) {
        mc->subtables[i] = malloc_or_die(N);
        memset(mc->subtables[i], 0, N);
    }

    return mc;
}


void markov_free(markov_t* mc)
{
    cond_dist16_free(&mc->catchall);
    size_t i;
    for (i = 0; i < NUM_SUBTABLES; ++i) {
        free(mc->subtables[i]);
    }
    free(mc);
}


/* Find a k-mer in the markov chain. */
static cell_t* markov_get(const markov_t* mc, kmer_t x)
{
    uint64_t h1, h0 = kmer_hash(x);
    uint16_t fp = fingerprint_hash(h0);

    uint64_t hs[NUM_SUBTABLES];

    size_t i, j;
    h1 = h0;
    for (i = 0; i < NUM_SUBTABLES; ++i) {
        h1 = kmer_hash_mix(h0, h1);
        hs[i] = h1 % mc->n;
        prefetch(&mc->subtables[i][NUM_CELLS_PER_BUCKET * hs[i]], 0, 0);
    }

    cell_t* c;
    for (i = 0; i < NUM_SUBTABLES; ++i) {
        c = &mc->subtables[i][NUM_CELLS_PER_BUCKET * hs[i]];
        for (j = 0; j < NUM_CELLS_PER_BUCKET; ++j) {
            if (c[j].fingerprint == 0) {
                c[j].fingerprint = fp;
                dist16_init(&c[j].dist);
                return &c[j];
            }
            else if (c[j].fingerprint == fp) {
                return &c[j];
            }
        }
    }

    return NULL;
}


/* Encode a k-mer with the given markov chain and update parameters. */
void markov_encode_and_update(markov_t* mc, ac_t* ac, kmer_t ctx, kmer_t x)
{
    kmer_t u = ctx & mc->xmask;
    cell_t* c = markov_get(mc, u);

    if (c == NULL) {
        cond_dist16_encode(ac, &mc->catchall, ctx & mc->xmask_catchall, x);
    }
    else {
        dist16_encode2(ac, &c->dist, x, 16);
    }
}


/* Decode a k-mer with the given markov chain and update parameters. */
kmer_t markov_decode_and_update(markov_t* mc, ac_t* ac, kmer_t ctx)
{
    kmer_t u = ctx & mc->xmask;
    cell_t* c = markov_get(mc, u);

    if (c == NULL) {
        return cond_dist16_decode(ac, &mc->catchall, ctx & mc->xmask_catchall);
    }
    else {
        return dist16_decode(ac, &c->dist);
    }
}


