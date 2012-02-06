
#include "qualenc.h"
#include "ac.h"
#include "dist.h"
#include "misc.h"
#include <stdlib.h>
#include <string.h>
#include <assert.h>


/* Every quality encoding scheme uses ASCII charocters in [33, 104] */
static const char   qual_last  = 40;
static const size_t qual_size  = 41;
static const size_t pos_bins   = 8;
static const size_t mu_bins    = 16;


struct qualenc_t_
{
    ac_t* ac;
    cond_dist41_t cs;
    dist41_t ms;
};


static void qualenc_init(qualenc_t* E, bool decode)
{
    cond_dist41_init(&E->cs, pos_bins * mu_bins * qual_size * qual_size, decode);
    dist41_init(&E->ms, decode);
}


qualenc_t* qualenc_alloc_encoder(quip_writer_t writer, void* writer_data)
{
    qualenc_t* E = malloc_or_die(sizeof(qualenc_t));
    E->ac = ac_alloc_encoder(writer, writer_data);

    qualenc_init(E, false);

    return E;
}


qualenc_t* qualenc_alloc_decoder(quip_reader_t reader, void* reader_data)
{
    qualenc_t* E = malloc_or_die(sizeof(qualenc_t));
    E->ac = ac_alloc_decoder(reader, reader_data);

    qualenc_init(E, true);

    return E;
}


void qualenc_free(qualenc_t* E)
{
    cond_dist41_free(&E->cs);
    dist41_free(&E->ms);
    ac_free(E->ac);
    free(E);
}


static inline char charmax2(char a, char b)
{
    return a > b ? a : b;
}



void qualenc_encode(qualenc_t* E, const seq_t* x)
{
    unsigned char q;   /* quality score */
    unsigned char q1;  /* preceding quality score */
    unsigned char q2;  /* maximum of three quality scores preceding q1 */

    char* qs = x->qual.s;
    size_t n = x->qual.n;

    size_t i;

    uint32_t mu = 0;
    for (i = 0; i < n; ++i) mu += qs[i];

    mu = (uint32_t) ((float) mu / (float) n);
    mu = (mu * mu_bins) / (qual_size + 1);

    dist41_encode(E->ac, &E->ms, mu);


    /* I am using the previous 4 positions to predict the next quality score.
     * These special cases here for positions 0 <= i <= 3 are a bit ugly but
     * speed the loop up considerably.
     */

    /* case: i = 0 */
    if (n >= 1) {
        cond_dist41_encode(E->ac, &E->cs,
                mu * qual_size * qual_size +
                qual_last * qual_size + qual_last, qs[0]);
    }

    /* case: i = 1 */
    if (n >= 2) {
        cond_dist41_encode(E->ac, &E->cs,
                mu * qual_size * qual_size +
                qual_last * qual_size + qs[0], qs[1]);
    }

    /* case: i = 2 */
    if (n >= 3) {
        cond_dist41_encode(E->ac, &E->cs,
                mu * qual_size * qual_size +
                qs[0] * qual_size + qs[1], qs[2]);
    }

    /* case: i >= 4 */
    for (i = 3; i < n; ++i) {
        q  = qs[i];
        q1 = qs[i - 1];
        q2 = charmax2(qs[i - 2], qs[i - 3]);

        cond_dist41_encode(E->ac, &E->cs,
                (i * pos_bins) / (n + 1) * mu_bins * qual_size * qual_size +
                mu * qual_size * qual_size +
                q2 * qual_size + q1, q);
    }
}


void qualenc_flush(qualenc_t* E)
{
    ac_flush_encoder(E->ac);
}


void qualenc_decode(qualenc_t* E, seq_t* seq, size_t n)
{
    str_t* qual = &seq->qual;
    while (n >= qual->size) fastq_expand_str(qual);
    qual->n = 0;
    char* qs = seq->qual.s;

    unsigned char q1;  /* preceding quality score */
    unsigned char q2;  /* maximum of three quality scores preceding q1 */

    uint32_t mu = dist41_decode(E->ac, &E->ms);

    size_t i;

    /* case: i = 0 */
    if (n >= 1) {
        qs[qual->n++] =
            cond_dist41_decode(E->ac, &E->cs,
                    mu * qual_size * qual_size +
                    qual_last * qual_size + qual_last);
    }

    /* case: i = 1 */
    if (n >= 2) {
        qs[qual->n++] =
            cond_dist41_decode(E->ac, &E->cs,
                    mu * qual_size * qual_size +
                    qual_last * qual_size + qs[0]);
    }

    /* case: i = 2 */
    if (n >= 3) {
        qs[qual->n++] =
            cond_dist41_decode(E->ac, &E->cs,
                    mu * qual_size * qual_size +
                    qs[0] * qual_size + qs[1]);
    }


    i = 3;
    while (qual->n < n) {
        q1 = qual->s[qual->n - 1];
        q2 = charmax2(
                qual->s[qual->n - 2],
                qual->s[qual->n - 3]);

        qs[qual->n++] =
            cond_dist41_decode(E->ac, &E->cs,
                    (i * pos_bins) / (n + 1) * mu_bins * qual_size * qual_size +
                    mu * qual_size * qual_size +
                    q2 * qual_size + q1);
        ++i;
    }

    qual->s[qual->n] = '\0';
}


void qualenc_reset_decoder(qualenc_t* E)
{
    ac_reset_decoder(E->ac);
}


