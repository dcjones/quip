
#include "qualenc.h"
#include "cumdist.h"
#include "misc.h"
#include <stdlib.h>
#include <string.h>
#include <assert.h>


/* Every quality encoding scheme uses ASCII charocters in [33, 104] */
const char   qual_first = 33;
const char   qual_last  = 104;
const size_t qual_size  = 72;

/* We expand arrays as we see longer reads. This is the expansion regime, which
 * is based on common read lengths*/
/*const size_t ms[] = {36, 50, 55, 75, 100, 150};*/
const size_t ms[] = {150};
const size_t ms_size = sizeof(ms);


/*
 */
typedef struct qualmodel_t_
{
    /* Score given position, nucleotide, and previous score */
    cumdist_t** cs;

    /* Maximum read length. */
    size_t m;
} qualmodel_t;



qualmodel_t* qualmodel_alloc()
{
    qualmodel_t* M = malloc_or_die(sizeof(qualmodel_t));
    M->m = ms[0];
    size_t i;

    M->cs = malloc_or_die(M->m * 25 * qual_size * sizeof(cumdist_t*));
    for (i = 0; i < M->m * 25 * qual_size; ++i) {
        M->cs[i] = cumdist_alloc(qual_size);
    }

    return M;
}



void qualmodel_free(qualmodel_t* M)
{
    size_t i;
    for (i = 0; i < M->m * 25 * qual_size; ++i) {
        cumdist_free(M->cs[i]);
    }
    free(M->cs);

    free(M);
}


unsigned char ascii_to_nucnum(char x)
{
    switch (x) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default:  return 4;
    }
}



struct qualenc_t_
{
    qualmodel_t* M;

    /* base and length of the current interval */
    uint32_t b, l;

    uint8_t* buf;
    size_t buflen;
    size_t bufpos; /* next vacant buffer position */

    quip_block_writer_t writer;
    void* writer_data;
};


qualenc_t* qualenc_alloc(quip_block_writer_t writer, void* writer_data)
{
    qualenc_t* E = malloc_or_die(sizeof(qualenc_t));
    E->M = qualmodel_alloc();

    E->writer      = writer;
    E->writer_data = writer_data;

    E->b = 0x00000000;
    E->l = 0xffffffff;

    E->bufpos = 0;
    E->buflen = 1024;
    E->buf    = malloc_or_die(E->buflen * sizeof(uint8_t));

    return E;
}


void qualenc_free(qualenc_t* E)
{
    qualmodel_free(E->M);
    free(E->buf);
    free(E);
}


static void qualenc_buf_append(qualenc_t* E, uint8_t c)
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


static void qualenc_encoder_renormalization(qualenc_t* E)
{
    uint32_t v;
    while (E->l < 0xffffff) {
        v = E->b >> 24;
        qualenc_buf_append(E, (uint8_t) v);
        E->l <<= 8;
        E->b <<= 8;
    }
}


static void qualenc_propogate_carry(qualenc_t* E)
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




static void qualenc_interval_update(qualenc_t* E, uint32_t p, uint32_t P)
{
    uint32_t y = ((uint64_t) E->l * (uint64_t) P) >> 32;
    uint32_t x = ((uint64_t) E->l * (uint64_t) (P - p)) >> 32;

    assert(y > x);

    uint32_t b0 = E->b;

    E->b += x;
    E->l = y - x;

    if (b0 > E->b) qualenc_propogate_carry(E);

    qualenc_encoder_renormalization(E);
}


/* encode a single quality/run-length pair */
static void qualenc_encode_qk(qualenc_t* E,
                              size_t        i,
                              unsigned char u,
                              unsigned char q0,
                              unsigned char q)
{
    uint64_t Z;
    uint32_t p, P;

    /* encode quality score */
    cumdist_t* cs = E->M->cs[i * (25 * qual_size) + u * qual_size + q0];
    Z = cumdist_Z(cs);
    p = ((uint64_t) cumdist_p(cs, q) << 32) / Z;
    P = ((uint64_t) cumdist_P(cs, q) << 32) / Z;

    if (P == 0) P = 0xffffffff;
    qualenc_interval_update(E, p, P);
}


void qualenc_encode(qualenc_t* E, const seq_t* x)
{
    if (x->qual.n == 0) return;

    if (x->qual.n > E->M->m) {
        // TODO: handle max read length increase
    }

    /* the number of preceeding quality scores to use to predict the next */
    static const size_t k = 3;

    unsigned char q0; /* preceeding quality score */
    unsigned char q;  /* quality score */
    unsigned char u;  /* nucleotide */

    size_t i, j;
    q0 = 0;
    for (i = 0; i < x->qual.n; ++i) {

        u = i < x->qual.n - 1 ? 5 * ascii_to_nucnum(x->seq.s[i + 1]) : 0;
        u += ascii_to_nucnum(x->seq.s[i]);

        q = x->qual.s[i] - qual_first;

        /* use the maximum of the last k quality scores to predict the next */
        q0 = qual_first;
        j = i < k ?  0 : i - k;
        for (; j < i; ++j) if (x->qual.s[j] > q0) q0 = x->qual.s[j];
        q0 -= qual_first;

        qualenc_encode_qk(E, i, u, q0, q);

        /* update model */
        cumdist_add(E->M->cs[i * (25 * qual_size) + u * qual_size + q0], q, 1);
    }
}


void qualenc_clear(qualenc_t* E)
{
    E->writer(E->writer_data, E->buf, E->bufpos);
    E->bufpos = 0;
}


