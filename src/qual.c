
#include "qual.h"
#include "cumdist.h"
#include "misc.h"
#include <stdlib.h>
#include <string.h>
#include <assert.h>


/* Every quality encoding scheme uses ASCII charocters in [33, 126] */
const char   qual_first = 33;
const char   qual_last  = 126;
const size_t qual_size  = 94;

/* We expand arrays as we see longer reads. This is the expansion regime, which
 * is based on common read lengths*/
/*const size_t ms[] = {36, 50, 55, 75, 100, 150};*/
const size_t ms[] = {100, 150};
const size_t ms_size = sizeof(ms);


/* A statistical model of quality scores. This is actually a pretty simple
 * model, so allow me to describe it here.
 *
 * Quality scores are essentially run-length encoded, but we assign a
 * propability to each run length and condition on position and score.
 * A probability is assigned to the quality score my conditioning on nucleotide
 * and position.
 *
 * Drawn as a graphical model, this looks something like:
 *
 *                 N    P
 *                 |   /|
 *                 |  / |
 *                 | /  |
 *                 |/   |
 *                 V    V
 *                 Q--->K
 *
 * Where N is the nucleotide, P is the position, Q is the quality score, and K
 * is the run-length.
 */
typedef struct qualmodel_t_
{
    /* Score given position and nucleotide. */
    cumdist_t** cs;

    /* Run-length given position and score. */
    cumdist_t** ck;

    /* Maximum read length. */
    size_t m;
} qualmodel_t;




qualmodel_t* qualmodel_alloc()
{
    qualmodel_t* M = malloc_or_die(sizeof(qualmodel_t));
    M->m = ms[0];
    size_t i, j;

    M->cs = malloc_or_die(M->m * 5 * sizeof(cumdist_t*));
    for (i = 0; i < M->m * 5; ++i) {
        M->cs[i] = cumdist_alloc(qual_size);
    }

    M->ck = malloc_or_die(M->m * qual_size * sizeof(cumdist_t*));
    for (i = 0; i < M->m; ++i) {
        for (j = 0; j < qual_size; ++j) {
            M->ck[i * qual_size + j] = cumdist_alloc(M->m - i);
        }
    }

    return M;
}



void qualmodel_free(qualmodel_t* M)
{
    size_t i;
    for (i = 0; i < M->m * 5; ++i) {
        cumdist_free(M->cs[i]);
    }
    free(M->cs);

    for (i = 0; i < M->m * qual_size; ++i) {
        cumdist_free(M->ck[i]);
    }
    free(M->ck);

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


void qualmodel_update(qualmodel_t* M, const seq_t* x)
{
    unsigned char q; // quality score
    unsigned char u; // nucleotide
    size_t ki, kj;   // position
    size_t l;        // run-length

    assert(x->qual.n <= M->m);

    u = ascii_to_nucnum(x->seq.s[0]);
    q = (unsigned char) (x->qual.s[0] - qual_first);
    cumdist_add(M->cs[0 * 5 + u], q, 1);

    ki = 0;
    kj = 1;
    while (kj < x->qual.n) {
        u = ascii_to_nucnum(x->seq.s[kj]);
        q = (unsigned char) (x->qual.s[kj] - qual_first);
        cumdist_add(M->cs[kj * 5 + u], q, 1);

        if (x->qual.s[kj] != x->qual.s[ki]) {
            l = kj - ki;
            q = (unsigned char) (x->qual.s[ki] - qual_first);
            cumdist_add(M->ck[ki * qual_size + q], l - 1, 1);

            ki = kj;
        }

        ++kj;
    }

    l = kj - ki;
    q = (unsigned char) (x->qual.s[ki] - qual_first);
    cumdist_add(M->ck[ki * qual_size + q], l - 1, 1);
}


struct qualenc_t_
{
    qualmodel_t* M;

    /* base and length of the current interval */
    uint32_t b, l;

    uint8_t* buf;
    size_t buflen;
    size_t bufpos; /* next vacant buffer position */

    qualenc_writer_t writer;
    void* writer_data;
};


qualenc_t* qualenc_alloc(qualenc_writer_t writer, void* writer_data)
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
    while (E->l < 0xffffff) {
        qualenc_buf_append(E, (uint8_t) (E->b >> 24));
        E->b <<= 8;
        E->l <<= 8;
        /*assert((uint64_t) E->b + (uint64_t) E->l < 0x100000000);*/
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

    /*assert((uint64_t) E->b + (uint64_t) x < 0x100000000);*/

    E->b += x;
    E->l = y - x;

    /*assert((uint64_t) E->b + (uint64_t) E->l < 0x100000000);*/

    if (b0 > E->b) qualenc_propogate_carry(E);

    qualenc_encoder_renormalization(E);
}


/* encode a single quality/run-length pair */
static void qualenc_encode_qk(qualenc_t* E,
                              size_t        pos,
                              size_t        len,
                              unsigned char nuc,
                              unsigned char qual)
{
    uint64_t Z;
    uint32_t p, P;

    /* encode quality score */
    cumdist_t* cs = E->M->cs[pos * 5 + nuc];
    Z = cumdist_Z(cs);
    p = ((uint64_t) cumdist_p(cs, qual) << 32) / Z;
    P = ((uint64_t) cumdist_P(cs, qual) << 32) / Z;
    if (P == 0) P = 0xffffffff;
    qualenc_interval_update(E, p, P);

    /* encode run-length */
    cumdist_t* ck = E->M->ck[pos * qual_size + qual];
    Z = cumdist_Z(ck);
    p = ((uint64_t) cumdist_p(ck, len) << 32) / Z;
    P = ((uint64_t) cumdist_P(ck, len) << 32) / Z;
    if (P == 0) P = 0xffffffff;
    qualenc_interval_update(E, p, P);
}


void qualenc_encode(qualenc_t* E, const seq_t* x)
{
    if (x->qual.n == 0) return;

    if (x->qual.n > E->M->m) {
        // TODO: handle max read length increase
    }

    unsigned char q; /* quality score */
    unsigned char u; /* nucleotide */
    size_t ki, kj;   /* position */
    size_t l;        /* run-length */

    ki = 0;
    kj = 1;
    do {
        if (kj >= x->qual.n || x->qual.s[kj] != x->qual.s[ki]) {
            l = kj - ki;
            u = ascii_to_nucnum(x->seq.s[ki]);
            q = (unsigned char) (x->qual.s[ki] - qual_first);
            qualenc_encode_qk(E, ki, l - 1, u, q);

            ki = kj;
        }

        ++kj;
    } while (kj < x->qual.n);

    qualmodel_update(E->M, x);
}


void qualenc_finish(qualenc_t* E)
{
    /* evict the buffer */
    E->writer(E->writer_data, E->buf, E->bufpos);
}


