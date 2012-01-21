
#include "ac.h"
#include "misc.h"
#include <stdlib.h>
#include <assert.h>
#include <string.h>

static const uint32_t min_length = 0x01000000U;
static const uint32_t max_length = 0xFFFFFFFFU;


struct ac_t_
{
    /* coder state */
    uint32_t b; /* base */
    uint32_t l; /* length */
    uint32_t v; /* value */

    uint8_t* buf;
    size_t buflen;
    size_t bufpos; /* next vacant buffer position */

    quip_writer_t writer;
    void* writer_data;

    quip_reader_t reader;
    void* reader_data;
};


ac_t* ac_alloc_encoder(quip_writer_t writer, void* writer_data)
{
    ac_t* ac = malloc_or_die(sizeof(ac_t));
    ac->b = 0;
    ac->l = max_length;

    ac->buflen = 4096;
    ac->buf = malloc_or_die(ac->buflen * sizeof(uint8_t));
    ac->bufpos = 0;

    ac->writer = writer;
    ac->writer_data = writer_data;

    ac->reader = NULL;
    ac->reader_data = NULL;

    return ac;
}


// TODO
/*ac_t* ac_alloc_decoder(quip_reader_t reader, void* reader_data);*/


void  ac_free(ac_t* ac)
{
    free(ac->buf);
    free(ac);
}



static void ac_append_byte(ac_t* E, uint8_t c)
{
    /* when the buffer is full, try to evict as much as possible, then shift
     * everything over */
    if (E->bufpos >= E->buflen) {
        size_t i = E->buflen - 1;
        while (i > 0) {
            if (E->buf[i] == 0xff && E->buf[i - 1] < 0xff) break;
            --i;
        }

        if (i > 1) {
            E->writer(E->writer_data, E->buf, i - 1);
            memmove(E->buf, E->buf + (i - 1), (E->buflen - (i - 1)) * sizeof(uint8_t));
            E->bufpos = E->buflen - (i - 1);
        }
        else if (i == 0) {
            E->writer(E->writer_data, E->buf, E->buflen - 1);
            E->buf[0] = E->buf[E->buflen - 1];
            E->bufpos = 1;
        }
    }

    /* if the buffer is _still_ full, we must expand it. */
    if (E->bufpos >= E->buflen) {
        E->buflen *= 2;
        E->buf = realloc_or_die(E->buf, E->buflen * sizeof(uint8_t));
     }

    assert(E->bufpos < E->buflen);

    E->buf[E->bufpos++] = c;
}


static void ac_propogate_carry(ac_t* ac)
{
    if (ac->bufpos == 0) return;

    size_t i;
    for (i = ac->bufpos - 1; i > 0; --i) {
        if (ac->buf[i] < 0xff) {
            ac->buf[i] += 1;
            break;
        }
        else ac->buf[i] = 0;
    }

    if (i == 0) {
        assert(ac->buf[0] < 0xff);
        ac->buf[0] += 1;
    }
}


static void ac_renormalize_encoder(ac_t* ac)
{
    while (ac->l < min_length) {
        ac_append_byte(ac, (uint8_t) (ac->b >> 24));
        ac->b <<= 8;
        ac->l <<= 8;
    }
}


void ac_encode(ac_t* ac, dist_t* D, symb_t x)
{
    uint32_t u, b0 = ac->b;

    if (x == D->n - 1) {
        u = D->ps[x] * (ac->l >> dist_length_shift);
        ac->b += u;
        ac->l -= u;
    }
    else {
        /* TODO: we might consider using SSE4.1 instructions to do these two
         * multiplications in one operation. */

        u = D->ps[x] * (ac->l >>= dist_length_shift);
        ac->b += u;
        ac->l = D->ps[x + 1] * ac->l - u;
    }

    if (b0 > ac->b)         ac_propogate_carry(ac);
    if (ac->l < min_length) ac_renormalize_encoder(ac);

    /* this is "dist_add(D, x, 1)", inlined for performance */
    D->cs[x]++;
    D->z++;
    if (--D->update_delay == 0) dist_update(D);
}


void ac_flush_encoder(ac_t* ac)
{
    uint32_t b0 = ac->b;

    if (ac->l > 2 * min_length) {
        ac->b += min_length;
        ac->l = min_length >> 1;
    }
    else {
        ac->b += min_length >> 1;
        ac->l = min_length >> 9;
    }

    if (b0 > ac->b) ac_propogate_carry(ac);
    ac_renormalize_encoder(ac);

}
