#include "bitset.h"
#include <assert.h>
#include <stdint.h>
#include <string.h>
#include <limits.h>

/* TODO:
 * All or most of this code assumes big-endian encoding of integers.
 * We may want to do something about that eventually.
 */


typedef uint_fast32_t bitset_int_t;

struct bitset_t_
{
    bitset_int_t* xs;
    size_t n;

};


static const size_t nbytes = sizeof(bitset_int_t);
static const size_t nbits  = CHAR_BIT * sizeof(bitset_int_t);
static const int onesbyte = (1 << (2 * CHAR_BIT)) - 1;


static inline size_t ints_needed(size_t bits)
{
    return bits / nbits + 1;
}


static inline size_t bytes_needed(size_t bits)
{
    return (bits / nbits + 1) * nbytes;
}


static inline bitset_int_t slot(size_t i) {
    return i / nbits;
}


static inline bitset_int_t mask(size_t i)
{
    return 1 << (i % nbits);
}


bitset_t* bitset_alloc(size_t n)
{
    bitset_t* bs = malloc(sizeof(bitset_t));
    assert(bs != NULL);

    bs->xs = malloc(bytes_needed(n));
    assert(bs->xs != NULL);

    bs->n  = n;

    bitset_alloff(bs);

    return bs;
}


void bitset_free(bitset_t* bs)
{
    if (bs != NULL) {
        free(bs->xs);
        free(bs);
    }
}


void bitset_copy(bitset_t* bs, const bitset_t* cs)
{
    if (ints_needed(bs->n) < ints_needed(cs->n)) {
        free(bs->xs);
        bs->xs = malloc(bytes_needed(cs->n));
        assert(bs->xs != NULL);
        bs->n = cs->n;
    }
    else if (ints_needed(bs->n) > ints_needed(cs->n)) {
        memset(bs->xs + ints_needed(cs->n),
               0, (ints_needed(cs->n) - ints_needed(bs->n)) * nbytes);
    }

    memcpy(bs->xs, cs->xs, bytes_needed(cs->n));
}


void bitset_flip(bitset_t* bs, size_t i)
{
    bs->xs[slot(i)] ^= mask(i);
}


void bitset_set(bitset_t* bs, size_t i, bool x)
{
    if (x) bitset_on(bs, i);
    else   bitset_off(bs, i);
}

bool bitset_get(const bitset_t* bs, size_t i)
{
    return (bool) (bs->xs[slot(i)] & mask(i));
}

void bitset_on(bitset_t* bs, size_t i)
{
    bs->xs[slot(i)] |= mask(i);
}

void bitset_off(bitset_t* bs, size_t i)
{
    bs->xs[slot(i)] &= ~mask(i);
}

void bitset_allon(bitset_t* bs)
{
    memset(bs->xs, onesbyte, bytes_needed(bs->n));
}

void bitset_alloff(bitset_t* bs)
{
    memset(bs->xs, 0, bytes_needed(bs->n));
}

void bitset_union(bitset_t* bs, const bitset_t* cs)
{
    size_t i;
    size_t m = bs->n < cs->n ? ints_needed(bs->n) : ints_needed(cs->n);

    for (i = 0; i < m; ++i) {
        bs->xs[i] |= cs->xs[i];
    }
}

void bitset_intersect(bitset_t* bs, const bitset_t* cs)
{
    size_t i;
    size_t m = bs->n < cs->n ? ints_needed(bs->n) : ints_needed(cs->n);
    for (i = 0; i < m; ++i) {
        bs->xs[i] &= cs->xs[i];
    }

    m = ints_needed(bs->n);
    for (i = 0; i < m; ++i) {
        bs->xs[i] = 0;
    }
}


size_t bitset_count_on(const bitset_t* bs)
{
    /* counting bits, the K&R way */
    size_t m = ints_needed(bs->n);
    size_t i;
    size_t c = 0;
    bitset_int_t x;
    for (i = 0; i < m; ++i) {
        x = bs->xs[i];
        while (x) {
            c++;
            x &= x - 1;
        }
    }

    return c;
}


void bitset_print(const bitset_t* bs, FILE* f)
{
    size_t i;
    for (i = 0; i < bs->n; ++i) {
        fputc(bitset_get(bs, i) ? 'X' : '.', f);
    }
}

