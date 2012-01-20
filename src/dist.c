
#include "dist.h"
#include "misc.h"
#include <string.h>
#include <stdbool.h>

struct dist_t_
{
    size_t n;

    /* symbol frequency */
    uint32_t* ps;

    /* cumulative symbol frequency */
    uint32_t* cs;

    /* symbol order, such that ps[ord[i]] is the freq. of symbol i */
    size_t* rev_ord;
    size_t* ord;

    /* xs[i] : number of occurance of symbol i */
    uint32_t* xs;

    /* total number of occurances */
    uint32_t z;

    /* has dist_add been called since the last call to dist_update */
    bool update_pending;
};


dist_t* dist_alloc(size_t n)
{
    dist_t* P = malloc_or_die(sizeof(dist_t));
    P->n = n;
    P->ps = malloc_or_die(n * sizeof(uint32_t));
    P->cs = malloc_or_die(n * sizeof(uint32_t));

    P->rev_ord = malloc_or_die(n * sizeof(size_t));
    P->ord     = malloc_or_die(n * sizeof(size_t));

    P->xs = malloc_or_die(n * sizeof(uint32_t));

    /* initialize to pseudocounts of 1 */
    size_t i;
    for (i = 0; i < n; ++i) P->xs[i] = 1;
    P->z = n;

    for (i = 0; i < n; ++i){
        P->ord[i] = P->rev_ord[i] = i;
    }

    P->update_pending = true;

    dist_update(P);

    return P;
}


void dist_free(dist_t* P)
{
    free(P->ps);
    free(P->cs);
    free(P->rev_ord);
    free(P->ord);
    free(P->xs);
    free(P);
}


uint32_t dist_P(const dist_t* P, size_t i)
{
    return P->ps[P->ord[i]];
}

uint32_t dist_C(const dist_t* P, size_t i)
{
    return P->cs[P->ord[i]];
}


void dist_add(dist_t* P, size_t i, uint32_t x)
{
    /* prevent precision loss in the arithmetic coder */
    if (P->z + x >= 0xffffff) return;

    P->xs[i] += x;
    P->z += x;
    P->update_pending = true;
}


static void sort_ord_insert(uint32_t* xs, size_t* ord, const size_t n)
{
    if (n < 2) return;

    uint32_t x;
    size_t o;
    int i, j;
    for (j = 1; j < (int) n; ++j) {
        x = xs[j];
        o = ord[j];
        for (i = j - 1; i >= 0 && x < xs[i]; --i) {
            xs[i + 1]  = xs[i];
            ord[i + 1] = ord[i];
        }

        xs[i + 1]  = x;
        ord[i + 1] = o;
    }
}

static void sort_ord_quick(uint32_t* xs, size_t* ord, const size_t n)
{
    if (n <= 10) {
        sort_ord_insert(xs, ord, n);
        return;
    }

    size_t pivot;
    /* choose pivot as the median of the first last and middle */
    if (xs[n/2] <= xs[n-1]) {
        if (xs[0] <= xs[n/2])      pivot = n / 2;
        else if (xs[0] <= xs[n-1]) pivot = 0;
        else                       pivot = n - 1;
    }
    else {
        if (xs[0] <= xs[n-1])      pivot = n - 1;
        else if (xs[0] <= xs[n/2]) pivot = 0;
        else                       pivot = n / 2;
    }

    uint32_t xs_tmp,  xs_piv  = xs[pivot];
    size_t   ord_tmp, ord_piv = ord[pivot];

    xs[pivot] = xs[n - 1];
    xs[n - 1] = xs_piv;

    ord[pivot] = ord[n - 1];
    ord[n - 1] = ord_piv;

    size_t i, j;
    for (i = 0, j = 0; j < n - 1; ++j) {
        if (xs[j] <= xs_piv) {
            xs_tmp = xs[i];
            xs[i]  = xs[j];
            xs[j]  = xs_tmp;

            ord_tmp = ord[i];
            ord[i]  = ord[j];
            ord[j]  = ord_tmp;

            ++i;
        }
    }

    xs[n - 1] = xs[i];
    xs[i] = xs_piv;

    ord[n - 1] = ord[i];
    ord[i] = ord_piv;

    sort_ord_quick(xs, ord, pivot);
    sort_ord_quick(xs + pivot + 1, ord + pivot + 1, n - pivot - 1);
}


/* sort an array and keep track of the order */
static void sort_ord_merge(uint32_t* xs, size_t* ord, const size_t n)
{
    if (n < 2) return;
    if (n <= 7) {
        sort_ord_insert(xs, ord, n);
        return;
    }

    sort_ord_merge(xs, ord, n / 2);
    sort_ord_merge(xs + n / 2, ord + n / 2, n - n / 2);

    /* merge */
    size_t i = 0;
    size_t j = n / 2;

    uint32_t xs_tmp;
    size_t ord_tmp;


    while (i < j && j < n) {
        if (xs[i] <= xs[j]) {
            ++i;
        }
        else {
            xs_tmp = xs[i];
            xs[i]  = xs[j];
            xs[j]  = xs_tmp;

            ord_tmp = ord[i];
            ord[i]  = ord[j];
            ord[j]  = ord_tmp;

            ++j;
        }
    }
}



void dist_update(dist_t* P)
{
    if (!P->update_pending) return;

    size_t i;
    for (i = 0; i < P->n; ++i) {
        P->ps[P->ord[i]] = P->xs[i];
    }

    /* sort symbols by frequency */

    /* If we have accumulated many samples, the order is unlikely to change, so
     * we use insertion sort, otherwise merge sort. */
    /*if (P->z < 10 * P->n && P->n > 7) {*/
        /*sort_ord_quick(P->ps, P->rev_ord, P->n);*/
    /*}*/
    /*else {*/
        sort_ord_insert(P->ps, P->rev_ord, P->n);
    /*}*/


    for (i = 0; i < P->n; ++i) {
        P->ord[P->rev_ord[i]] = i;
    }

    /* normalize */
    uint64_t Z = P->z;
    for (i = 0; i < P->n; ++i) {
        P->ps[i] = (uint32_t) (((uint64_t) P->ps[i] << 32) / Z);
    }

    /* accumulate */
    P->cs[0] = P->ps[0];
    for (i = 1; i < P->n - 1; ++i) {
        P->cs[i] = P->cs[i - 1] + P->ps[i];
    }
    P->cs[P->n - 1] = 0xffffffff;

    P->update_pending = false;
}


size_t dist_find(dist_t* P, uint32_t p)
{
    size_t a = 0;
    size_t b = P->n - 1;
    size_t c;

    while (a + 1 < b) {
        c = (a + b) / 2;

        if (P->cs[c] < p) a = c;
        else              b = c;
    }

    return P->rev_ord[b];
}

