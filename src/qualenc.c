
#include "qualenc.h"
#include "qualenc_prior.h"
#include "ac.h"
#include "dist.h"
#include "misc.h"
#include <stdlib.h>
#include <string.h>
#include <assert.h>

/* The rate at which the quality markov chain is updated. */
static const size_t qual_update_rate = 4;

/* Every quality encoding scheme uses ASCII charocters in [33, 104] */
static const char   qual_last  = 40;
static const size_t qual_size  = 41;
static const size_t pos_bins   = 4;
static const size_t q_bins     = 16;
static const size_t delta_bins = 8;
static const int    delta_max  = 50;

/* Map quality scores to a smaller alphabet size. */
static const uint8_t q_bin_map[41] =
  {  15, 15, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14,
     13, 13, 13, 13, 13, 13, 13, 13, 12, 12, 12, 12, 12, 12, 11, 10,
      9,  8,  7,  6,  5,  4,  3,  2,  1 };


/* Map running deltas to a smaller alphabet size. */
static const uint8_t delta_bin_map[50] =
  {   0,  1,  2,  2,  2,  3,  3,  3,  3,  3,
      4,  4,  4,  4,  4,  4,  4,  4,  5,  5,
      5,  5,  5,  5,  5,  5,  5,  5,  6,  6,
      6,  6,  6,  6,  6,  6,  6,  6,  6,  6,
      6,  6,  6,  6,  6,  6,  6,  6,  6,  7 };


struct qualenc_t_
{
    ac_t* ac;
    cond_dist41_t cs;
};


/* Compute indexes into the quality conditional distribution */
#if 0
#define cs_index(bin, delta, q3, q2, q1) \
 (q1 + qual_size * (\
    q2 + q_bins * (\
      q3 + q_bins * (\
        delta_bin_map[delta] + delta_bins * (\
          bin)))))
#endif

/* This version depends on the specific values defined above, but is somewhat faster. */
#define cs_index(bin, delta, q3, q2, q1) \
    ((q1) + qual_size * (((((((bin) << 3) | delta_bin_map[delta]) << 4) | (q3)) << 4) | (q2)))


static void qualenc_setprior(qualenc_t* E)
{
    size_t q0, q1;
    for (q1 = 0; q1 < qual_size; ++q1) {
        for (q0 = 0; q0 < qual_size; ++q0) {
            E->cs.xss[q1].xs[q0].count = qualenc_prior[q1 * qual_size + q0];
        }
        dist41_update(&E->cs.xss[q1]);
    }

    size_t i;
    for (i = 1; i < q_bins * q_bins * delta_bins * pos_bins; ++i) {
        memcpy(&E->cs.xss[i * qual_size], &E->cs.xss[0],
                qual_size * sizeof(dist41_t));
    }
}


static void qualenc_init(qualenc_t* E)
{
    cond_dist41_init(&E->cs, delta_bins * pos_bins * q_bins * q_bins * qual_size);
    qualenc_setprior(E);

    cond_dist41_set_update_rate(&E->cs, qual_update_rate);
}


qualenc_t* qualenc_alloc_encoder(quip_writer_t writer, void* writer_data)
{
    qualenc_t* E = malloc_or_die(sizeof(qualenc_t));
    E->ac = ac_alloc_encoder(writer, writer_data);

    qualenc_init(E);

    return E;
}


qualenc_t* qualenc_alloc_decoder(quip_reader_t reader, void* reader_data)
{
    qualenc_t* E = malloc_or_die(sizeof(qualenc_t));
    E->ac = ac_alloc_decoder(reader, reader_data);

    qualenc_init(E);

    return E;
}


void qualenc_free(qualenc_t* E)
{
    cond_dist41_free(&E->cs);
    ac_free(E->ac);
    free(E);
}

static inline char charmin2(char a, char b)
{
    return a < b ? a : b;
}


static inline char charmax2(char a, char b)
{
    return a > b ? a : b;
}

static inline int intmin2(int a, int b)
{
    return a < b ? a : b;
}


void qualenc_encode(qualenc_t* E, const seq_t* x)
{
    prefetch(x->qual.s);

    union {
        uint64_t ui64;
        uint8_t  ui8[4];
    } qprev;

    qprev.ui64 = 0;

    int delta = 0; 
    int qdiff;

    char* qs = x->qual.s;
    size_t n = x->qual.n;

    /* this is: ceil(n / pos_bins) */
    size_t pos_bin_size = (n + pos_bins - 1) / pos_bins;
    size_t i;

    for (i = 0; i < n; ++i) {
        cond_dist41_encode(E->ac, &E->cs,
                cs_index(i / pos_bin_size, delta,
                         charmin2(qprev.ui8[3], qprev.ui8[2]),
                         qprev.ui8[1], qprev.ui8[0]), qs[i]);

        qdiff = (int) qprev.ui8[0] - (int) qs[i];

        qprev.ui64 <<= 8;
        qprev.ui8[1] = q_bin_map[qprev.ui8[1]];
        qprev.ui8[0] = qs[i];

        if (qdiff != 0) {
            delta += 1;
            if (delta >= delta_max) {
                delta = delta_max - 1;
                ++i;
                break;
            }
        }
    }

    /* We enter this loop when delta hits its maximum. It is the same as the
     * above loop but saves a few cycles by not dealing with delta any further.
     * */
    for (; i < n; ++i) {
        cond_dist41_encode(E->ac, &E->cs,
                cs_index(i / pos_bin_size, delta,
                         charmin2(qprev.ui8[3], qprev.ui8[2]),
                         qprev.ui8[1], qprev.ui8[0]), qs[i]);

        qprev.ui64 <<= 8;
        qprev.ui8[1] = q_bin_map[qprev.ui8[1]];
        qprev.ui8[0] = qs[i];
    }
}

size_t qualenc_finish(qualenc_t* E)
{
    return ac_finish_encoder(E->ac);
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

    union {
        uint64_t ui64;
        uint8_t  ui8[4];
    } qprev;

    qprev.ui64 = 0;

    int delta = 0; 
    int qdiff;

    /* this is: ceil(n / pos_bins) */
    size_t pos_bin_size = (n + pos_bins - 1) / pos_bins;
    size_t i;

    for (i = 0; i < n; ++i) {
        qs[i] = cond_dist41_decode(E->ac, &E->cs,
                    cs_index(i / pos_bin_size, delta,
                             charmin2(qprev.ui8[3], qprev.ui8[2]),
                             qprev.ui8[1], qprev.ui8[0]));

        qdiff = (int) qprev.ui8[0] - (int) qs[i];

        qprev.ui64 <<= 8;
        qprev.ui8[1] = q_bin_map[qprev.ui8[1]];
        qprev.ui8[0] = qs[i];

        if (qdiff != 0) {
            delta += 1;
            if (delta >= delta_max) {
                delta = delta_max - 1;
                ++i;
                break;
            }
        }
    }

    /* Faster loop once delta hits its maximum. */
    for (; i < n; ++i) {
        qs[i] = cond_dist41_decode(E->ac, &E->cs,
                    cs_index(i / pos_bin_size, delta,
                             charmin2(qprev.ui8[3], qprev.ui8[2]),
                             qprev.ui8[1], qprev.ui8[0]));

        qprev.ui64 <<= 8;
        qprev.ui8[1] = q_bin_map[qprev.ui8[1]];
        qprev.ui8[0] = qs[i];
    }

    qual->n = n;
    qual->s[qual->n] = '\0';
}


void qualenc_start_decoder(qualenc_t* E)
{
    ac_start_decoder(E->ac);
}


void qualenc_reset_decoder(qualenc_t* E)
{
    ac_reset_decoder(E->ac);
}


