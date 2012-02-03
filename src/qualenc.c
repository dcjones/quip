
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
static const size_t pos_bins   = 8;

struct qualenc_t_
{
    ac_t* ac;
    cond_dist72_t cs;
};


static void qualenc_init(qualenc_t* E, bool decode)
{
    cond_dist72_init(&E->cs, pos_bins * qual_size * qual_size, decode);
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
    cond_dist72_free(&E->cs);
    ac_free(E->ac);
    free(E);
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

    char* qs = x->qual.s;
    size_t n = x->qual.n;

    size_t i;

    /* I am using the previous 4 positions to predict the next quality score.
     * These special cases here for positions 0 <= i <= 3 are a bit ugly but
     * speed the loop up considerably.
     */

    /* case: i = 0 */
    if (n >= 1) {
        cond_dist72_encode(E->ac, &E->cs,
                71 * qual_size + 71, qs[0]);
    }

    /* case: i = 1 */
    if (n >= 2) {
        cond_dist72_encode(E->ac, &E->cs,
                71 * qual_size + qs[0], qs[1]);
    }

    /* case: i = 2 */
    if (n >= 3) {
        cond_dist72_encode(E->ac, &E->cs,
                qs[0] * qual_size + qs[1], qs[2]);
    }

    /* case: i = 3 */
    if (n >= 4) {
        cond_dist72_encode(E->ac, &E->cs,
                charmax2(qs[0], qs[1]) * qual_size + qs[2], qs[3]);
    }

    /* case: i >= 4 */
    q2 = 0;
    for (i = 4; i < n; ++i) {
        q  = qs[i];
        q1 = qs[i - 1];
        q2 = charmax3(
                qs[i - 2],
                qs[i - 3],
                qs[i - 4]);

        cond_dist72_encode(E->ac, &E->cs,
                (i * pos_bins) / (n + 1) * qual_size * qual_size + q2 * qual_size + q1, q);
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

    size_t i;

    /* case: i = 0 */
    if (n >= 1) {
        qs[qual->n++] =
            cond_dist72_decode(E->ac, &E->cs,
                    71 * qual_size + 71);
    }

    /* case: i = 1 */
    if (n >= 2) {
        qs[qual->n++] =
            cond_dist72_decode(E->ac, &E->cs,
                    71 * qual_size + qs[0]);
    }

    /* case: i = 2 */
    if (n >= 3) {
        qs[qual->n++] =
            cond_dist72_decode(E->ac, &E->cs,
                    qs[0] * qual_size + qs[1]);
    }

    /* case: i = 3 */
    if (n >= 4) {
        qs[qual->n++] =
            cond_dist72_decode(E->ac, &E->cs,
                    charmax2(qs[0], qs[1]) * qual_size + qs[2]);
    }


    i = 4;
    while (qual->n < n) {
        q1 = qual->s[qual->n - 1];
        q2 = charmax3(
                qual->s[qual->n - 2],
                qual->s[qual->n - 3],
                qual->s[qual->n - 4]);

        qs[qual->n++] =
            cond_dist72_decode(E->ac, &E->cs,
                    (i * pos_bins) / (n + 1) * qual_size * qual_size + q2 * qual_size + q1);
        ++i;
    }

    qual->s[qual->n] = '\0';
}


void qualenc_reset_decoder(qualenc_t* E)
{
    ac_reset_decoder(E->ac);
}


