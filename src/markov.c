
#include <math.h>
#include <string.h>

#include "dist.h"
#include "markov.h"
#include "misc.h"

static const size_t max_count = 1 << 15;


/* compute ceil(a/b) without casting to floats */
static size_t size_t_ceil(size_t a, size_t b)
{
    return (a + b - 1) / b;
}


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
    for (i = 0; i < 16; ++i) {
        z += D->xs[i].count;
    }

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
        c += D->xs[i].count;
    }

    D->update_delay = 4;
}


void maj_dist16_init(maj_dist16_t* D)
{
    D->majority = 0;
    size_t i;
    for (i = 0; i < 16; ++i) D->xs[i].count = 1;
    maj_dist16_update(D);
    D->update_delay = 2;
}


void maj_dist16_encode(ac_t* ac, maj_dist16_t* D, kmer_t x)
{
    prefetch(D, 1, 0);

    uint32_t b0 = ac->b;

    uint32_t u;
    if (x == 15) {
        u = (uint32_t) D->xs[x].freq * (ac->l >> dist_length_shift);
        ac->b += u;
        ac->l -= u;
    }
    else {
        u = (uint32_t) D->xs[x].freq * (ac->l >>= dist_length_shift);
        ac->b += u;
        ac->l = (uint32_t) D->xs[x + 1].freq * ac->l - u;
    }

    if (b0 > ac->b)         ac_propogate_carry(ac);
    if (ac->l < min_length) ac_renormalize_encoder(ac);

    D->xs[x].count += 100;
    if (D->xs[x].count > D->xs[D->majority].count) {
        D->majority = x;
    }

    if(!--D->update_delay) maj_dist16_update(D);
}


kmer_t maj_dist16_decode(ac_t* ac, maj_dist16_t* D)
{
    // TODO: No! We can't do this since we need maj_dist16_update to be called.
    // */
    return dist16_decode(ac, (dist16_t*) D);
}


/* Number of sub-tables used by the dlCFB. */
#define NUM_SUBTABLES 4

/* Number of cells per bucket. */
#define NUM_CELLS_PER_BUCKET 8


/* Each cell holds a 17-bit fingerprint (hash) of the k-mer as well as a
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
    mc->n = size_t_ceil(n, NUM_CELLS_PER_BUCKET * NUM_SUBTABLES);
    mc->k = k;
    mc->xmask = kmer_mask(k);
    mc->k_catchall = k_catchall;
    mc->xmask_catchall = kmer_mask(k_catchall);

    cond_dist16_init(&mc->catchall, 1 << (2 * k_catchall));
    cond_dist16_set_update_rate(&mc->catchall, 2);

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


#if 0
/* Return true iff the k-mer is present. */
static bool markov_has(const markov_t* mc, kmer_t x)
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

    size_t bucket_sizes[NUM_SUBTABLES];
    cell_t* cells[NUM_SUBTABLES];

    cell_t* c;
    for (i = 0; i < NUM_SUBTABLES; ++i) {
        c = &mc->subtables[i][NUM_CELLS_PER_BUCKET * hs[i]];
        for (j = 0; j < NUM_CELLS_PER_BUCKET && c[j].fingerprint != 0; ++j) {
            if (c[j].fingerprint == fp) {
                return true;
            }
        }

        cells[i] = &c[j];
        bucket_sizes[i] = j;
    }

    return false;
}
#endif


/* Find a k-mer in the markov chain. */
static cell_t* markov_get(const markov_t* mc, kmer_t x, bool* new_cell)
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

    size_t bucket_sizes[NUM_SUBTABLES];
    cell_t* cells[NUM_SUBTABLES];

    cell_t* c;
    for (i = 0; i < NUM_SUBTABLES; ++i) {
        c = &mc->subtables[i][NUM_CELLS_PER_BUCKET * hs[i]];
        for (j = 0; j < NUM_CELLS_PER_BUCKET && c[j].fingerprint != 0; ++j) {
            if (c[j].fingerprint == fp) {
                if (new_cell) *new_cell = false;
                return &c[j];
            }
        }

        cells[i] = &c[j];
        bucket_sizes[i] = j;
    }

    size_t i_min = NUM_SUBTABLES;
    size_t min_bucket_size = NUM_CELLS_PER_BUCKET;
    for (i = 0; i < NUM_SUBTABLES && min_bucket_size > 0; ++i) {
        if (bucket_sizes[i] < min_bucket_size) {
            i_min = i;
            min_bucket_size = bucket_sizes[i];
        }
    }

    if (i_min < NUM_SUBTABLES) {
        cells[i_min]->fingerprint = fp;
        maj_dist16_init(&cells[i_min]->dist);
        if (new_cell) *new_cell = true;
        return cells[i_min];
    }
    else return NULL;
}


/* Return true if path coercion seems like a good move. */
static bool path_coercion_heuristic(const maj_dist16_t* d, kmer_t x)
{
    if ((int) d->xs[d->majority].freq > 4 * (int) d->xs[x].freq) return true;
    else return false;
}


kmer_t markov_encode_and_update(markov_t* mc, ac_t* ac, size_t i, kmer_t ctx, kmer_t x)
{
    static unsigned long N = 0;
    static unsigned long new_kmer_count = 0;
    bool new_kmer;

    kmer_t u = ctx & mc->xmask;
    cell_t* c;
    if (i < mc->k || (c = markov_get(mc, u, &new_kmer)) == NULL) {
        cond_dist16_encode(ac, &mc->catchall, ctx & mc->xmask_catchall, x);
        return x;
    }

    ++N;
    if (new_kmer) ++new_kmer_count;
    if (N % 1000000 == 0) {
        fprintf(stderr, "new_kmer_count = %0.1f%%\n",
                100.0 * (double) new_kmer_count / (double) N);
        N = 0;
        new_kmer_count = 0;
    }

    maj_dist16_encode(ac, &c->dist, x);
    if (path_coercion_heuristic(&c->dist, x)) {
        return c->dist.majority;
    }
    else return x;
}


/* TODO: rewrite the decode function */

/* Decode a k-mer with the given markov chain and update parameters. */
kmer_t markov_decode_and_update(markov_t* mc, ac_t* ac, kmer_t ctx)
{
    kmer_t u = ctx & mc->xmask;
    cell_t* c = markov_get(mc, u, NULL);

    if (c == NULL) {
        return cond_dist16_decode(ac, &mc->catchall, ctx & mc->xmask_catchall);
    }
    else {
        return maj_dist16_decode(ac, &c->dist);
    }
}


