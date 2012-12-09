
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

    D->xs[x].count += 16;
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

    cond_dist16_t d_delta;
    dist2_t d_type;

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
    cond_dist16_init(&mc->d_delta, 16);
    dist2_init(&mc->d_type);

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


/* A heuristic guess at how many fewer bits are needed to encode the natural
 * path versus the majority path */
static int natural_path_advantage(const markov_t* mc, const maj_dist16_t* d,
                                  kmer_t ctx, kmer_t x)
{
    UNUSED(mc);
    UNUSED(ctx);
    uint16_t m = d->majority;

#if 0
    bool hasx = markov_has(mc, mc->xmask & ((ctx << 4) | x));
    bool hasm = markov_has(mc, mc->xmask & ((ctx << 4) | m));

    if ((!hasx && hasm) && d->xs[m].count > d->xs[x].count) return -1;
    else return 1;

    /*if ((long) d->xs[m].count > (long) d->xs[x].count) return -1;*/
    /*return 1;*/
#endif

#if 0
    uint32_t fx = x < 15 ? d->xs[x + 1].freq - d->xs[x].freq :
                          0x8000 - d->xs[x].freq;
    uint32_t fm = m < 15 ? d->xs[m + 1].freq - d->xs[m].freq :
                          0x8000 - d->xs[m].freq;

    if (fm > 2 * fx) return -1;
    else return 1;
#endif

    if ((int) d->xs[m].count > 2 * (int) d->xs[x].count) return -1;
    else return 1;
}



kmer_t markov_encode_and_update(markov_t* mc, ac_t* ac, size_t i, kmer_t ctx, kmer_t x)
{
    static uint64_t branch_a_count = 0;
    static uint64_t branch_b_count = 0;
    static uint64_t branch_c_count = 0;
    static uint64_t branch_d_count = 0;
    static uint64_t N = 0;
    static uint64_t new_cells = 0;

    ++N;

    if (i < mc->k) {
        cond_dist16_encode(ac, &mc->catchall, ctx & mc->xmask_catchall, x);
        return x;
    }

    kmer_t u = ctx & mc->xmask;
    bool new_cell;
    cell_t* c = markov_get(mc, u, &new_cell);
    if (new_cell) {
        ++new_cells;
    }

    if ((branch_a_count + branch_b_count + branch_c_count + branch_d_count) % 1000000 == 0) {
        fprintf(stderr, "%lu, %lu, %lu, %lu | %lu / %lu\n",
                (unsigned long) branch_a_count,
                (unsigned long) branch_b_count,
                (unsigned long) branch_c_count,
                (unsigned long) branch_d_count,
                (unsigned long) new_cells,
                (unsigned long) N);
        N = 0;
        new_cells = 0;
    }

    if (c == NULL) {
        ++branch_a_count;
        cond_dist16_encode(ac, &mc->catchall, ctx & mc->xmask_catchall, x);
        return x;
    }
    else {

        /* x in on the consensus path */
        if (x == c->dist.majority) {
            /*c->dist.xs[c->dist.majority].count += 16;*/
            /*if(!--c->dist.update_delay) maj_dist16_update(&c->dist);*/
            maj_dist16_encode(ac, &c->dist, x);
            ++branch_b_count;
            return x;
        }
        else if (natural_path_advantage(mc, &c->dist, ctx, x) >= 0) {
            /*c->dist.xs[x].count += 16;*/
            /*if(!--c->dist.update_delay) maj_dist16_update(&c->dist);*/

            /*cond_dist16_encode(ac, &mc->d_delta, x, c->dist.majority);*/

            maj_dist16_encode(ac, &c->dist, x);
            ++branch_c_count;
            return x;
        }
        /* encode a coercion. */
        else {
            maj_dist16_encode(ac, &c->dist, x);
            ++branch_d_count;
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
    cell_t* c = markov_get(mc, u, NULL);

    if (c == NULL) {
        return cond_dist16_decode(ac, &mc->catchall, ctx & mc->xmask_catchall);
    }
    else {
        return maj_dist16_decode(ac, &c->dist);
    }
}


