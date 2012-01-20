
#include "cumdist.h"
#include "misc.h"
#include <string.h>
#include <assert.h>

struct cumdist_t_
{
    /* subtree frequencies */
    uint32_t* fs;

    /* left subtree frequencies */
    uint32_t* ls;

    /* number of elements */
    size_t n;
};



static inline size_t parent_idx (size_t i) { return (i - 1) / 2; }
static inline size_t left_idx   (size_t i) { return 2 * i + 1; }
static inline size_t right_idx  (size_t i) { return 2 * i + 2; }

/* index of the i'th leaf node */
static inline size_t leaf_idx (size_t i, size_t n)
{
    return 2 * n - 2 - i;
}


/* If all (or many) of the leaves change at once, call this function to
 * efficiently update the internal nodes of the heap. */
static void cumdist_update_tree(cumdist_t* C)
{
    size_t i, j;

    for (i = leaf_idx(C->n - 1, C->n) - 1; 1; --i) {
        C->fs[i] = 0;

        j = left_idx(i);
        if (j < 2 * C->n - 1) {
            C->ls[i] = C->fs[j];
            C->fs[i] += C->fs[j];
        }

        j = right_idx(i);
        if (j < 2 * C->n - 1) {
            C->fs[i] += C->fs[j];
        }

        if (i == 0) break;
    }
}


cumdist_t* cumdist_alloc(size_t n)
{
    cumdist_t* C = malloc_or_die(sizeof(cumdist_t));
    C->n = n;

    C->fs = malloc_or_die((2 * n - 1) * sizeof(uint32_t));
    memset(C->fs, 0, (2 * n - 1) * sizeof(uint32_t));

    C->ls = malloc_or_die((n - 1) * sizeof(uint32_t));
    memset(C->ls, 0, (n - 1) * sizeof(uint32_t));

    /* initialize every element to a pseudocount of one */
    size_t i;
    for (i = 0; i < n; ++i) {
        C->fs[leaf_idx(i, C->n)] = 1;
    }

    cumdist_update_tree(C);

    return C;
}


void cumdist_free(cumdist_t* C)
{
    free(C->fs);
    free(C->ls);
    free(C);
}


size_t cumdist_n(const cumdist_t* C) { return C->n; }


uint32_t cumdist_p(const cumdist_t* C, size_t i)
{
    /* ith leaf index */
    i = 2 * C->n - 2 - i;
    return C->fs[i];
}


uint32_t cumdist_p_norm(const cumdist_t* C, size_t i)
{
    /* ith leaf index */
    i = 2 * C->n - 2 - i;
    return ((uint64_t) C->fs[i] << 32) / C->fs[0];
}


uint32_t cumdist_c(const cumdist_t* C, size_t i)
{
    /* ith leaf index */
    i = 2 * C->n - 2 - i;

    uint32_t c = C->fs[i];
    while (i > 0) {
        if (i % 2 == 0) c += C->ls[parent_idx(i)];
        i = parent_idx(i);
    }

    return c;
}


uint32_t cumdist_c_norm(const cumdist_t* C, size_t i)
{
    /* ith leaf index */
    i = 2 * C->n - 2 - i;

    uint32_t c = C->fs[i];
    while (i > 0) {
        if (i % 2 == 0) c += C->ls[parent_idx(i)];
        i = parent_idx(i);
    }

    c = ((uint64_t) c << 32) / C->fs[0];
    if (c == 0) c = 0xffffffff;

    return c;
}


uint32_t cumdist_z(const cumdist_t* C)
{
    return C->fs[0];
}

void cumdist_check(cumdist_t* C)
{
    size_t N = 2 * C->n - 1;
    uint32_t c;
    size_t i;
    for (i = 0; i < N - C->n; ++i) {
        c = 0;
        if (left_idx(i) < N) {
            c += C->fs[left_idx(i)];
            assert(C->fs[left_idx(i)] == C->ls[i]);
        }

        if (right_idx(i) < N) {
            c += C->fs[right_idx(i)];
        }

        assert(C->fs[i] == c);
    }
}


void cumdist_add(cumdist_t* C, size_t i, uint32_t x)
{
    /* prevent precision loss in the arithmetic coder */
    if (C->fs[0] >= 0xffffff) return;

    /* ith leaf index */
    i = 2 * C->n - 2 - i;

    C->fs[i] += x;

    while (i > 0) {
        if (i % 2 == 1) C->ls[parent_idx(i)] += x;
        i = parent_idx(i);
        C->fs[i] += x;
    }
}


void cumdist_expand(cumdist_t* C, size_t new_n)
{
    assert(new_n >= C->n);

    cumdist_t* D = cumdist_alloc(new_n);
    size_t i;
    for (i = 0; i < C->n; ++i) {
        cumdist_add(D, i, cumdist_p(C, i) - 1);
    }

    uint32_t* tmp;

    tmp = C->fs;
    C->fs = D->fs;
    D->fs = tmp;

    tmp = C->ls;
    C->ls = D->ls;
    D->ls = tmp;

    D->n = C->n;
    C->n = new_n;

    cumdist_free(D);
}


