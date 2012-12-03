
#include <math.h>
#include <string.h>

#include "dist.h"
#include "markov.h"
#include "misc.h"

#if 0
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
#endif

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
    for (i = 0; i < 16; ++i) {
        if (D->xs[i].count > D->xs[D->majority].count) D->majority = i;
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

    D->xs[x].count += 8;

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

    /* Encode the observed dinucleotide, conditioned on the majority
     * dinucleotide, when doing path coercion */
    cond_dist16_t d_delta;

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
    cond_dist16_init(&mc->d_delta, 16);

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


/* A heuristic guess at how many fewer bits are needed to encode the natural
 * path versus the majority path */
static int natural_path_advantage(const markov_t* mc, const maj_dist16_t* d,
                                  kmer_t x)
{
    UNUSED(mc);
    uint16_t m = d->majority;
    uint32_t fx = x < 15 ? d->xs[x + 1].freq - d->xs[x].freq :
                          0x8000 - d->xs[x].freq;
    uint32_t fm = m < 15 ? d->xs[m + 1].freq - d->xs[m].freq :
                          0x8000 - d->xs[m].freq;

    if (fm > 10 * fx) return -1;
    else return 1;
}



kmer_t markov_encode_and_update(markov_t* mc, ac_t* ac, kmer_t ctx, kmer_t x)
{
    static uint64_t branch_a_count = 0;
    static uint64_t branch_b_count = 0;
    static uint64_t branch_c_count = 0;

    kmer_t u = ctx & mc->xmask;
    cell_t* c = markov_get(mc, u);

    if (c == NULL) {
        ++branch_a_count;
        cond_dist16_encode(ac, &mc->catchall, ctx & mc->xmask_catchall, x);
        return x;
    }
    else {
        /* x in on the consensus path */
        if (x == c->dist.majority || natural_path_advantage(mc, &c->dist, x) >= 0) {
            ++branch_b_count;
            dist2_encode(ac, &mc->d_natural_coerced_path, MARKOV_TYPE_NATURAL);
            maj_dist16_encode(ac, &c->dist, x);
            return x;
        }
        /* encode a coercion. */
        else {
            ++branch_c_count;
            dist2_encode(ac, &mc->d_natural_coerced_path, MARKOV_TYPE_COERCION);
            cond_dist16_encode(ac, &mc->d_delta, c->dist.majority, x);
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


