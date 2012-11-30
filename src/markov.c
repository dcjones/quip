
#include <math.h>
#include <string.h>

#include "dist.h"
#include "markov.h"
#include "misc.h"

/* Number of leading zeros in a uint16_t. */
static int nlz(uint16_t x)
{
    int n = 16;
    uint16_t y;
    y = x >> 8; if (y != 0) {n -= 8; x = y;}
    y = x >> 4; if (y != 0) {n -= 4; x = y;}
    y = x >> 2; if (y != 0) {n -= 2; x = y;}
    y = x >> 1; if (y != 0) return n - 2;
    return n - x;
}


/* (Very) approximate -log2(p) for fixed-precision probability stored in a uint16_t. */
static int inlog2(uint16_t p)
{
    return 0x1000000 - nlz(p);
}


static const size_t max_count = 1 << 15;

/* This structure is essentially dist16_t, but with an added field that tracks
 * the highest probability element, which is necessary to implement heuristics
 * for finding good paths in the de Bruijn graph. */
typedef struct maj_dist16_t_
{
    uint16_t update_delay;

    struct {
        uint16_t count;
        uint16_t freq;
    } xs[16];

    /* Keep this here so we can cast this to dist16_t */
    uint16_t majority;
} maj_dist16_t;


void maj_dist16_update(maj_dist16_t* D)
{
    size_t i;

    /* rescale when we have exceeded the maximum count */
    uint32_t z = 0;
    for (i = 0; i < 16; ++i) z += D->xs[i].count;
    if (z > max_count) {
        z = 0;
        for (i = 0; i < 16; ++i) {
            D->xs[i].count /= 2;
            if (D->xs[i].count == 0) D->xs[i].count = 1;
            z += D->xs[i].count;
        }
    }


    /* update frequencies */
    const uint32_t scale = 0x80000000U / z;
    const uint32_t shift = 31 - dist_length_shift;
    uint32_t c = 0;

    for (i = 0; i < 16; ++i) {
        D->xs[i].freq = (uint16_t) ((scale * c) >> shift);
        if (D->xs[i].freq > D->xs[D->majority].freq) D->majority = i;
        c += D->xs[i].count;
    }

    D->update_delay = 16;
}


void maj_dist16_init(maj_dist16_t* D)
{
    D->update_delay = 16;
    D->majority = 0;
    size_t i;
    for (i = 0; i < 16; ++i) D->xs[i].count = 1;
    maj_dist16_update(D);
}


void maj_dist16_encode(ac_t* ac, maj_dist16_t* D, kmer_t x)
{
    dist16_encode(ac, (dist16_t*) D, x);
}


kmer_t maj_dist16_decode(ac_t* ac, maj_dist16_t* D)
{
    return dist16_decode(ac, (dist16_t*) D);
}


/* Number of sub-tables used by the dlCFB. */
#define NUM_SUBTABLES 4

/* Number of cells per bucket. */
#define NUM_CELLS_PER_BUCKET 8


/* Each cell holds a 16-bit fingerprint (hash) of the k-mer as well as a
 * conditional distribution over nucleotides following that k-mer. */
typedef struct cell_t_
{
    uint16_t fingerprint;
    maj_dist16_t dist;
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


/* Are we taking the "natural" edge dictated by the sequence, or coercing it
 * down the majority path? */
enum {
    MARKOV_TYPE_NATURAL,
    MARKOV_TYPE_COERCION
};


/* The idea here is to use the dlCBF to map a k-mer to a distribution. */
struct markov_t_
{
    cell_t* subtables[NUM_SUBTABLES];

    /* If a k-mer could not be found or inserted in the dlcbf, encode against
     * a lower-order dense markov chain. */
    cond_dist16_t catchall;

    /* Encode a bit that determines whether wether or not we are making an
     * alteration to the sequence being encoded in order to follow the majority
     * path. */
    dist2_t d_natural_coerced_path;

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
    dist2_init(&mc->d_natural_coerced_path);

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
                maj_dist16_init(&c[j].dist);
                return &c[j];
            }
            else if (c[j].fingerprint == fp) {
                return &c[j];
            }
        }
    }

    return NULL;
}


kmer_t markov_encode_and_update(markov_t* mc, ac_t* ac, kmer_t ctx, kmer_t x)
{
    kmer_t u = ctx & mc->xmask;
    cell_t* c = markov_get(mc, u);

    if (c == NULL) {
        cond_dist16_encode(ac, &mc->catchall, ctx & mc->xmask_catchall, x);
        return x;
    }
    else {
        /* x in on the consensus path */
        if (x == c->dist.majority ||
            mc->k * inlog2(c->dist.xs[x].freq) <=
            mc->k * inlog2(c->dist.xs[c->dist.majority].freq) + 4) {
            dist2_encode(ac, &mc->d_natural_coerced_path, MARKOV_TYPE_NATURAL);
            maj_dist16_encode(ac, &c->dist, x);
            return x;
        }
        /* Encode a modification. */
        else {
            dist2_encode(ac, &mc->d_natural_coerced_path, MARKOV_TYPE_COERCION);
            /* We are always coercing to the majority path, so no more
             * information is needed. */
            /* TODO: we may want to update the count on the majority path */
            return c->dist.majority;
        }
    }
}


#if 0
/* Encode a k-mer with the given markov chain and update parameters. */
void markov_encode_and_update(markov_t* mc, ac_t* ac, kmer_t ctx, kmer_t x)
{
    kmer_t u = ctx & mc->xmask;
    cell_t* c = markov_get(mc, u);

    if (c == NULL) {
        cond_dist16_encode(ac, &mc->catchall, ctx & mc->xmask_catchall, x);
    }
    else {
        maj_dist16_encode(ac, &c->dist, x);
    }
}
#endif


/* TODO: rewrite the decode function */

/* Decode a k-mer with the given markov chain and update parameters. */
kmer_t markov_decode_and_update(markov_t* mc, ac_t* ac, kmer_t ctx)
{
    kmer_t u = ctx & mc->xmask;
    cell_t* c = markov_get(mc, u);

    if (c == NULL) {
        return cond_dist16_decode(ac, &mc->catchall, ctx & mc->xmask_catchall);
    }
    else {
        return maj_dist16_decode(ac, &c->dist);
    }
}


