
#include "qualenc.h"
#include "ac.h"
#include "dist.h"
#include "misc.h"
#include <stdlib.h>
#include <string.h>
#include <assert.h>

/* The rate at which the quality markov chain is updated. */
static const size_t qual_update_rate = 6;

/* These constants should not be tweaked willy-nilly. The tables
   below and the cs_index macro depend on them being what they are. */
static const size_t qual_size  = 64;
static const size_t pos_bins   = 4;
static const size_t q_bins1    = 42;
static const size_t q_bins2    = 16;
static const size_t delta_bins = 8;
static const int    delta_max  = 50;

/* Map a quality score alphabet of size 64 to an alphabet of size of 42. */
static const uint8_t q_bin_map1[64] =
  {  1,  2,  3,  4,  5,  6,  7,  8,
     9, 10, 11, 12, 13, 14, 15, 16,
    17, 18, 19, 20, 21, 22, 23, 24,
    25, 26, 27, 28, 29, 30, 31, 32,
    33, 34, 35, 36, 38, 39, 40, 41,
    41, 41, 41, 41, 41, 41, 41, 41,
    41, 41, 41, 41, 41, 41, 41, 41,
    41, 41, 41, 41, 41, 41, 41, 41 };

/* Map a quality score alphabet of size 42 to an alphabet of size of 16. */
static const uint8_t q_bin_map2[42] =
  {  1,  1,  2,  2,  2,  2,  2,  2,
     2,  2,  3,  3,  3,  3,  3,  3,
     3,  3,  3,  3,  3,  3,  3,  3,
     3,  4,  4,  4,  4,  4,  5,  6,
     7,  8,  9, 10, 11, 12, 13, 14,
    15, 15 };

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
    cond_dist64_t cs;
    uint8_t base_qual;
};


/* Compute an index into the conditional distribution over quality scores. */
#define cs_index(pos_bin, delta, q3, q2, q1) \
    ((((((((pos_bin << 3) | \
     delta_bin_map[delta]) << 4) | \
      (q3)) << 4) | \
       (q2)) * 42) + \
        (q1))


static void qualenc_init(qualenc_t* E)
{
    cond_dist64_init(&E->cs, delta_bins * pos_bins * q_bins2 * q_bins2 * q_bins1);
    cond_dist64_set_update_rate(&E->cs, qual_update_rate);
    E->base_qual = '!';
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
    cond_dist64_free(&E->cs);
    ac_free(E->ac);
    free(E);
}


void qualenc_set_base_qual(qualenc_t* E, char base_qual)
{
    E->base_qual = base_qual;
}


static inline uint8_t bytemax2(uint8_t a, uint8_t b)
{
    return a > b ? a : b;
}

static inline int intmin2(int a, int b)
{
    return a < b ? a : b;
}


void qualenc_encode(qualenc_t* E, const short_read_t* x)
{
    union {
        uint64_t ui64;
        uint8_t  ui8[4];
    } qprev;

    qprev.ui64 = 0;

    int delta = 0; 
    int qdiff;
    uint8_t q;

    uint8_t* qs = x->qual.s;
    size_t n = x->qual.n;

    /* this is: ceil(n / pos_bins) */
    size_t pos_bin_size = (n + pos_bins - 1) / pos_bins;
    size_t i;

    for (i = 0; i < n; ++i) {
        q = qs[i] - E->base_qual;
        cond_dist64_encode(E->ac, &E->cs,
                cs_index(i / pos_bin_size, delta,
                         bytemax2(qprev.ui8[3], qprev.ui8[2]),
                         qprev.ui8[1], qprev.ui8[0]), q);

        qdiff = (int) qprev.ui8[0] - (int) q;

        qprev.ui64 <<= 8;
        qprev.ui8[1] = q_bin_map2[qprev.ui8[1]];
        qprev.ui8[0] = q_bin_map1[q];

        if (qdiff < -1 || qdiff >= 1) {
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
        q = qs[i] - E->base_qual;
        cond_dist64_encode(E->ac, &E->cs,
                cs_index(i / pos_bin_size, delta,
                         bytemax2(qprev.ui8[3], qprev.ui8[2]),
                         qprev.ui8[1], qprev.ui8[0]), q);

        qprev.ui64 <<= 8;
        qprev.ui8[1] = q_bin_map2[qprev.ui8[1]];
        qprev.ui8[0] = q_bin_map1[q];
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


void qualenc_decode(qualenc_t* E, short_read_t* seq, size_t n)
{
    str_t* qual = &seq->qual;
    str_reserve(qual, n + 1);
    qual->n = 0;
    uint8_t* qs = seq->qual.s;

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
        qs[i] = cond_dist64_decode(E->ac, &E->cs,
                    cs_index(i / pos_bin_size, delta,
                             bytemax2(qprev.ui8[3], qprev.ui8[2]),
                             qprev.ui8[1], qprev.ui8[0]));

        qdiff = (int) qprev.ui8[0] - (int) qs[i];

        qprev.ui64 <<= 8;
        qprev.ui8[1] = q_bin_map2[qprev.ui8[1]];
        qprev.ui8[0] = q_bin_map1[qs[i]];
        qs[i] += E->base_qual;

        if (qdiff < -1 || qdiff >= 1) {
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
        qs[i] = cond_dist64_decode(E->ac, &E->cs,
                    cs_index(i / pos_bin_size, delta,
                             bytemax2(qprev.ui8[3], qprev.ui8[2]),
                             qprev.ui8[1], qprev.ui8[0]));

        qprev.ui64 <<= 8;
        qprev.ui8[1] = q_bin_map2[qprev.ui8[1]];
        qprev.ui8[0] = q_bin_map1[qs[i]];
        qs[i] += E->base_qual;
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


