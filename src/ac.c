
#include "ac.h"
#include "misc.h"
#include <stdlib.h>
#include <assert.h>
#include <string.h>


struct ac_t_
{
    uint64_t b, l;

    uint8_t* buf;
    size_t buflen;
    size_t bufpos; /* next vacant buffer position */

    quip_block_writer_t writer;
    void* writer_data;
};


static void ac_buf_append(ac_t* E, uint8_t c)
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


static void ac_renormalize(ac_t* E)
{
    uint32_t v;
    while (E->l < 0xffffff) {
        v = E->b >> 24;
        ac_buf_append(E, (uint8_t) v);
        E->l <<= 8;
        E->b <<= 8;
    }
}


static void ac_propogate_carry(ac_t* E)
{
    size_t i;
    for (i = E->bufpos - 1; i > 0; --i) {
        if (E->buf[i] < 0xff) {
            E->buf[i] += 1;
            break;
        }
        else E->buf[0] = 0;
    }

    if (i == 0) {
        assert(E->buf[0] < 0xff);
        E->buf[0] += 1;
    }
}



ac_t* ac_alloc(quip_block_writer_t writer, void* writer_data)
{
    ac_t* E = malloc_or_die(sizeof(ac_t));

    E->writer = writer;
    E->writer_data = writer_data;

    E->b = 0x00000000;
    E->l = 0xffffffff;

    E->bufpos = 0;
    E->buflen = 1024;
    E->buf = malloc_or_die(E->buflen * sizeof(uint8_t));

    return E;
}


void ac_free(ac_t* E)
{
    free(E->buf);
    free(E);
}



void ac_update(ac_t* E, uint64_t p, uint64_t c)
{
    uint64_t y = (E->l * c) >> 32;
    uint64_t x = (E->l * (c - p)) >> 32;

    assert(y > x);

    uint64_t b0 = E->b;

    E->b += x; // XXX mod 2^32 ??
    E->l = y - x;

    if (b0 > E->b) ac_propogate_carry(E);

    ac_renormalize(E);
}


void ac_flush(ac_t* E)
{
    uint32_t b0 = E->b;

    E->b = E->b + 0x800000; /* b = b + D^(P-1) / 2 */
    E->l = 0xffff;          /* l = D^(P - 2) - 1 */

    if (b0 > E->b) ac_propogate_carry(E);

    ac_renormalize(E);


    /* now, reset for continued use */

    E->b = 0x00000000;
    E->l = 0xffffffff;
    E->bufpos = 0;

}


struct dec_t_
{
    quip_reader_t reader;
    void* reader_data;

    /* interval length */
    uint64_t l;

    /* intput buffer */
    uint8_t* buf;
    size_t bufsize;

    /* the next byte in the input buffer */
    uint8_t next;

    /* available bytes */
    size_t avail;
};


static void dec_refill(dec_t* D)
{
    if (D->avail > 0 && D->next != D->buf) {
        memmove(D->buf, D->next, D->avail);
        D->next = D->buf;
    }

    if (D->avail < D->bufsize) {
        D->avail += D->reader(D->reader_data, D->buf + D->avail, D->bufsize - D->avail);
    }
}


dec_t* dec_alloc(quip_reader_t reader, void* reader_data)
{
    dec_t* D = malloc_or_die(sizeof(dec_t));

    D->reader      = reader;
    D->reader_data = reader_data;
    D->l = 0xffffffff;

    D->bufsize = 4096;
    D->buf = malloc_or_die(D->bufsize * sizeof(uint8_t));
    D->avail = 0;
    D->next = D->buf;
}


void dec_free(dec_t* D)
{
    free(D->buf);
    free(D);
}




