
#include "qualenc.h"
#include "ac.h"
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
        case 'A': case 'a': case 'U': case 'u': return 1;
        case 'C': case 'c': return 2;
        case 'G': case 'g': return 3;
        case 'T': case 't': return 4;
        default:  return 0;
    }
}



struct qualenc_t_
{
    qualmodel_t* M;
    ac_t* ac;
};


qualenc_t* qualenc_alloc(quip_block_writer_t writer, void* writer_data)
{
    qualenc_t* E = malloc_or_die(sizeof(qualenc_t));
    E->M  = qualmodel_alloc();
    E->ac = ac_alloc(writer, writer_data);

    return E;
}


void qualenc_free(qualenc_t* E)
{
    qualmodel_free(E->M);
    ac_free(E->ac);
    free(E);
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
    ac_update(E->ac, p, P);
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


void qualenc_flush(qualenc_t* E)
{
    // TODO
}


