
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

static const size_t qual_size_sq = 72 * 72;

/* Dependence on position within a read is modeled by using a fixed number of
 * bins. */
static const size_t read_pos_bins = 10;


struct qualenc_t_
{
    ac_t* ac;
    dist_t** cs;
};


static void qualenc_init(qualenc_t* E, bool decode)
{
    E->cs = dist_alloc_array(read_pos_bins * qual_size * qual_size, qual_size, decode);
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
    dist_free_array(E->cs);
    ac_free(E->ac);
    free(E);
}


static void qualenc_encode_pos(qualenc_t* E,
                              size_t        i,
                              unsigned char q2,
                              unsigned char q1,
                              unsigned char q)
{
    size_t idx = i  * qual_size_sq +
                 q1 * qual_size +
                 q2;

    dist_t* cs = E->cs[idx];
    ac_encode(E->ac, cs, q);
}


static unsigned char qualenc_decode_pos(
            qualenc_t* E,
            size_t        i,
            unsigned char q2,
            unsigned char q1)
{
    size_t idx = i  * qual_size_sq +
                 q1 * qual_size +
                 q2;

    dist_t* cs = E->cs[idx];
    return ac_decode(E->ac, cs);
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

    size_t bin_size = (n / read_pos_bins) + 1;
    size_t i;


    /* I am using the previous 4 positions to predict the next quality score.
     * These special cases here for positions 0 <= i <= 3 are a bit ugly but
     * speed the loop up considerably.
     */

    /* case: i = 0 */
    if (n >= 1) {
        qualenc_encode_pos(E, 0 / bin_size,
                0, 0, qs[0]);
    }

    /* case: i = 1 */
    if (n >= 2) {
        qualenc_encode_pos(E, 1 / bin_size,
                0, qs[0], qs[1]);
    }

    /* case: i = 2 */
    if (n >= 3) {
        qualenc_encode_pos(E, 2 / bin_size,
                qs[0], qs[1], qs[2]);
    }

    /* case: i = 3 */
    if (n >= 4) {
        qualenc_encode_pos(E, 3 / bin_size,
                charmax2(qs[0], qs[1]),
                qs[2], qs[3]);
    }

    /* case: i >= 4 */
    for (i = 4; i < n; ++i) {

        q  = qs[i];
        q1 = qs[i - 1];
        q2 = charmax3(
                qs[i - 2],
                qs[i - 3],
                qs[i - 4]);

        qualenc_encode_pos(E, i / bin_size, q2, q1, q);
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

    unsigned char q1;  /* preceding quality score */
    unsigned char q2;  /* maximum of three quality scores preceding q1 */

    size_t bin_size = (n / read_pos_bins) + 1;

    /* case: i = 0 */
    if (n >= 1) {
        qual->s[qual->n++] = qualenc_decode_pos(E, 0 / bin_size, 0, 0);
    }

    /* case: i = 1 */
    if (n >= 2) {
        qual->s[qual->n++] = qualenc_decode_pos(
                E, 1 / bin_size, 0, qual->s[0]);
    }

    /* case: i = 2 */
    if (n >= 3) {
        qual->s[qual->n++] = qualenc_decode_pos(
                E, 2 / bin_size, qual->s[0], qual->s[1]);
    }

    /* case: i = 3 */
    if (n >= 4) {
        qual->s[qual->n++] = qualenc_decode_pos(
                E, 3 / bin_size, charmax2(qual->s[0], qual->s[1]), qual->s[2]);
    }


    while (qual->n < n) {
        q1 = qual->s[qual->n - 1];
        q2 = charmax3(
                qual->s[qual->n - 2],
                qual->s[qual->n - 3],
                qual->s[qual->n - 4]);

        qual->s[qual->n] = qualenc_decode_pos(E, qual->n / bin_size, q2, q1);
        ++qual->n;
    }

    qual->s[qual->n] = '\0';
}


void qualenc_reset_decoder(qualenc_t* E)
{
    ac_reset_decoder(E->ac);
}


