
#include "qualenc.h"
#include "ac.h"
#include "dist.h"
#include "misc.h"
#include <stdlib.h>
#include <string.h>
#include <assert.h>


/* Every quality encoding scheme uses ASCII charocters in [33, 104] */
static const char   qual_first = 33;
static const char   qual_last  = 104;
static const size_t qual_size  = 72;

/* We expand arrays as we see longer reads. This is the expansion regime, which
 * is based on common read lengths*/
/*const size_t ms[] = {36, 50, 55, 75, 100, 150};*/
static const size_t ms[] = {150};
static const size_t ms_size = sizeof(ms);

/* Dependence on position within a read is modeled by using a fixed number of
 * bins. */
static const size_t read_pos_bins = 10;

/* How many dist_add calls are made before all the distributions are updated. */
static const size_t dist_update_delay = 25000;


/*
 */
typedef struct qualmodel_t_
{
    /* Score given position, nucleotide, and previous score */
    dist_t** cs;

    /* Maximum read length. */
    size_t m;
} qualmodel_t;



qualmodel_t* qualmodel_alloc()
{
    qualmodel_t* M = malloc_or_die(sizeof(qualmodel_t));
    M->m = ms[0] / read_pos_bins;
    size_t i;

    M->cs = malloc_or_die(M->m * qual_size * qual_size * sizeof(dist_t*));
    for (i = 0; i < M->m * qual_size * qual_size; ++i) {
        M->cs[i] = dist_alloc(qual_size);
    }

    return M;
}



void qualmodel_free(qualmodel_t* M)
{
    size_t i;
    for (i = 0; i < M->m * qual_size * qual_size; ++i) {
        dist_free(M->cs[i]);
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
    size_t update_delay;
    ac_t* ac;
};


qualenc_t* qualenc_alloc(quip_block_writer_t writer, void* writer_data)
{
    qualenc_t* E = malloc_or_die(sizeof(qualenc_t));
    E->M  = qualmodel_alloc();
    E->ac = ac_alloc(writer, writer_data);
    E->update_delay = dist_update_delay;

    return E;
}


void qualenc_free(qualenc_t* E)
{
    qualmodel_free(E->M);
    ac_free(E->ac);
    free(E);
}


static void qualenc_update_dist(qualenc_t* E)
{
    size_t i;
    for (i = 0; i < E->M->m * qual_size * qual_size; ++i) {
        dist_update(E->M->cs[i]);
    }
}

/* encode a single quality/run-length pair */
static void qualenc_encode_qk(qualenc_t* E,
                              size_t        i,
                              unsigned char q2,
                              unsigned char q1,
                              unsigned char q)
{
    uint32_t p, c;

    size_t idx = (i / read_pos_bins) * qual_size * qual_size +
                 q1 * qual_size +
                 q2;

    /* encode quality score */
    dist_t* cs = E->M->cs[idx];

    p = dist_P(cs, q);
    c = dist_C(cs, q);

    ac_update(E->ac, p, c);

    dist_add(cs, q, 1);
}


static inline char charmax2(char a, char b)
{
    return a > b ? a : b;
}

static inline char charmax3(char a, char b, char c)
{
    char d = a > b ? a : b;
    return c > d ? c : d;
}


void qualenc_encode(qualenc_t* E, const seq_t* x)
{
    unsigned char q;   /* quality score */
    unsigned char q1;  /* preceding quality score */
    unsigned char q2;  /* maximum of three quality scores preceding q1 */

    size_t i;

    /* I am using the previous 4 positions to predict the next quality score.
     * These special cases here for positions 0 <= i <= 3 are a bit ugly but
     * speed the loop up considerably.
     */

    /* case: i = 0 */
    if (x->qual.n >= 1) {
        qualenc_encode_qk(E, 0, 0, 0, x->qual.s[0]);
    }

    /* case: i = 1 */
    if (x->qual.n >= 2) {
        qualenc_encode_qk(E, 1, 0, x->qual.s[0], x->qual.s[1]);
    }

    /* case: i = 2 */
    if (x->qual.n >= 3) {
        qualenc_encode_qk(E, 2, x->qual.s[0], x->qual.s[1], x->qual.s[2]);
    }

    /* case: i = 3 */
    if (x->qual.n >= 4) {
        qualenc_encode_qk(E, 3,
                charmax2(x->qual.s[0], x->qual.s[1]),
                x->qual.s[2], x->qual.s[3]);
    }

    /* case: i >= 4 */
    for (i = 4; i < x->qual.n; ++i) {

        q  = x->qual.s[i];
        q1 = x->qual.s[i - 1];
        q2 = charmax3(
                x->qual.s[i - 2],
                x->qual.s[i - 3],
                x->qual.s[i - 4]);


        qualenc_encode_qk(E, i, q2, q1, q);
    }

    if (--E->update_delay == 0) {
        qualenc_update_dist(E);
        E->update_delay = dist_update_delay;
    }
}


void qualenc_flush(qualenc_t* E)
{
    ac_flush(E->ac);
}


