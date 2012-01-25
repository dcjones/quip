
#include "ac.h"
#include "misc.h"
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <stdbool.h>

/* allowable values of the interval length before renormalization */
static const uint32_t min_length = 0x01000000U;
static const uint32_t max_length = 0xFFFFFFFFU;


struct ac_t_
{
    /* coder state */
    uint32_t b; /* base */
    uint32_t l; /* length */
    uint32_t v; /* value */

    /* input or output buffer (depending on whether we are decoding or encoding,
     * respectively). */
    uint8_t* buf;

    /* size allocated to buf */
    size_t buflen;

    /* index next vacant buffer position */
    size_t bufpos;

    /* available indput (for decoding) */
    size_t bufavail;

    /* callback function for encoder output */
    quip_writer_t writer;
    void* writer_data;

    /* callback function for decoder input */
    quip_reader_t reader;
    void* reader_data;

    /* initial state: used when decoding */
    bool init_state;
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


ac_t* ac_alloc_decoder(quip_reader_t reader, void* reader_data)
{
    ac_t* ac = malloc_or_die(sizeof(ac_t));
    ac->l = max_length;

    ac->buflen = 4096;
    ac->buf = malloc_or_die(ac->buflen * sizeof(uint8_t));
    ac->bufpos = 0;

    ac->writer = NULL;
    ac->writer_data = NULL;

    ac->reader = reader;
    ac->reader_data = reader_data;

    ac->init_state = true;

    return ac;
}


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


static uint8_t ac_get_byte(ac_t* ac)
{
    if (ac->bufpos < ac->bufavail) {
        return ac->buf[ac->bufpos++];
    }
    else {
        ac->bufavail = ac->reader(ac->reader_data, ac->buf, ac->buflen);
        ac->bufpos = 0;
        if (ac->bufpos < ac->bufavail) return ac->buf[ac->bufpos++];
        else return 0;
    }
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


static void ac_renormalize_decoder(ac_t* ac)
{
    while (ac->l < min_length) {
        ac->v = (ac->v << 8) | (uint32_t) ac_get_byte(ac);
        ac->l <<= 8;
    }
}


void ac_encode(ac_t* ac, dist_t* D, symb_t x)
{
    uint32_t u, b0 = ac->b;

    u = D->ps[x] * (ac->l >>= dist_length_shift);
    ac->b += u;
    ac->l = D->ps[x + 1] * ac->l - u;

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

    ac->writer(ac->writer_data, ac->buf, ac->bufpos);

    /* reset encoder */
    ac->b = 0;
    ac->l = max_length;
    ac->bufpos = 0;
}


symb_t ac_decode(ac_t* ac, dist_t* D)
{
    if (ac->init_state) {
        ac->bufavail = ac->reader(ac->reader_data, ac->buf, ac->buflen);

        if (ac->bufavail < 4) {
            fprintf(stderr, "Malformed compressed data encountered.");
            exit(EXIT_FAILURE);
        }

        ac->v = ((uint32_t) ac->buf[0] << 24) | ((uint32_t) ac->buf[1] << 16) |
                ((uint32_t) ac->buf[2] << 8)  | ((uint32_t) ac->buf[3]);

        ac->l = max_length;

        ac->bufpos = 4;
        ac->init_state = false;
    }

    symb_t s, n, m;
    uint32_t z, x, y = ac->l;

    if (D->dec) {
        uint32_t dv = ac->v / (ac->l >>= dist_length_shift);
        uint32_t t = dv >> D->dec_shift;

        s = D->dec[t];
        n = D->dec[t + 1] + 1;

        /* binary search in [s, n] */
        while (n > s + 1) {
            m = (s + n) >> 1;
            if (D->ps[m] > dv) n = m;
            else               s = m;
        }

        x = D->ps[s] * ac->l;
        y = D->ps[s + 1] * ac->l;
    }
    else {
        /* decode using binary search */
        x = s = 0;
        ac->l >>= dist_length_shift;
        n = D->n;
        m = n >> 1;

        do {
            z = ac->l * D->ps[m];
            if (z > ac->v) {
                n = m;
                y = z;
            }
            else {
                s = m;
                x = z;
            }

            m = (s + n) >> 1;

        } while(m != s);
    }

    ac->v -= x;
    ac->l = y - x;

    if (ac->l < min_length) ac_renormalize_decoder(ac);

    /* this is "dist_add(D, x, 1)", inlined for performance */
    D->cs[s]++;
    D->z++;
    if (--D->update_delay == 0) dist_update(D);

    return s;
}


void ac_reset_decoder(ac_t* ac)
{
    ac->bufpos = 0;
    ac->bufavail = 0;
    ac->init_state = true;
}


