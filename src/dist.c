
#include "dist.h"
#include "misc.h"
#include <string.h>

/* The distribution is updated every "initial_update_factor * n" new
 * observations. */
static const size_t update_delay_factor = 1;
const size_t dist_length_shift = 15;
static const size_t max_count = 1 << 15;


dist_t* dist_alloc(size_t n, bool decode)
{
    return decode ? dist_alloc_decode(n) : dist_alloc_encode(n);
}


dist_t* dist_alloc_encode(size_t n)
{
    dist_t* D = malloc_or_die(sizeof(dist_t));
    D->n = n;
    D->last_symbol = n - 1;
    D->ps = malloc_or_die(n * sizeof(uint32_t));
    D->cs = malloc_or_die(n * sizeof(uint32_t));
    D->dec = NULL;
    D->update_delay = D->n * update_delay_factor;

    /* initialize to pseudocounts of 1 */
    size_t i;
    for (i = 0; i < n; ++i) D->cs[i] = 1;
    D->z = n;

    dist_update(D);

    return D;
}


dist_t* dist_alloc_decode(size_t n)
{
    dist_t* D = malloc_or_die(sizeof(dist_t));
    D->n = n;
    D->last_symbol = n - 1;
    D->ps = malloc_or_die(n * sizeof(uint32_t));
    D->cs = malloc_or_die(n * sizeof(uint32_t));
    D->update_delay = D->n * update_delay_factor;

    if (n > 16) {
        size_t dec_bits = 3;
        while (n > (1U << (dec_bits + 2))) ++dec_bits;
        D->dec_size  = (1 << dec_bits) + 4;
        D->dec_shift = dist_length_shift - dec_bits;
        D->dec = malloc_or_die((D->dec_size + 6) * sizeof(uint32_t));
    }
    else {
        D->dec_shift = 0;
        D->dec_size  = 0;
        D->dec = NULL;
    }

    /* initialize to pseudocounts of 1 */
    size_t i;
    for (i = 0; i < n; ++i) D->cs[i] = 1;
    D->z = n;

    dist_update(D);

    return D;
}


void dist_free(dist_t* D)
{
    if (D == NULL) return;
    free(D->ps);
    free(D->cs);
    free(D->dec);
    free(D);
}


size_t dist_n(const dist_t* D)
{
    return D->n;
}


uint32_t dist_P(const dist_t* D, symb_t x)
{
    return D->ps[x];
}


void dist_add(dist_t* D, symb_t x, uint32_t k)
{
    D->cs[x] += k;
    D->z += k;

    if (--D->update_delay == 0) {
        dist_update(D);
    }
}

void dist_update(dist_t* D)
{
    size_t i;

    /* rescale when we have exceeded the maximum count */
    if (D->z > max_count) {
        D->z = 0;
        for (i = 0; i < D->n; ++i) {
            D->cs[i] = D->cs[i] / 2 + 1;
            D->z += D->cs[i];
        }
    }

    /* update frequencies */
    const uint32_t scale = 0x80000000U / D->z;
    const uint32_t shift = 31 - dist_length_shift;
    uint32_t j, w, c = 0;

    for (i = 0; i < D->n; ++i) {
        D->ps[i] = (scale * c) >> shift;
        c += D->cs[i];
    }

    /* update decoder table */
    if (D->dec) {
        for (i = 0, j = 0; i < D->n; ++i) {
            w = D->ps[i] >> D->dec_shift;
            while (j < w) D->dec[++j] = i - 1;
        }
        D->dec[0] = 0;
        while (j < D->dec_size) D->dec[++j] = D->n - 1;
    }


    D->update_delay = D->n * update_delay_factor;
}

