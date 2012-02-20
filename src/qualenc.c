
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
static const size_t pos_bins   = 6;
static const size_t q2_bins    = 12;
static const size_t q3_bins    = 12;
static const size_t delta_bins = 8;
static const int    delta_max  = 41;


struct qualenc_t_
{
    ac_t* ac;
    cond_dist41_t cs;
};


#define cs_index(n, i, delta, q3, q2, q1) \
    (q1 + qual_size * (\
            (qual_size - q2 >= q2_bins ? q2_bins - 1 : qual_size - q2) + q2_bins * (\
                (qual_size - q3 >= q3_bins ? q3_bins - 1 : qual_size - q3) + q3_bins * (\
                    (delta * delta_bins / delta_max) + delta_bins * (\
                        (i * pos_bins) / n)))))



static void qualenc_init(qualenc_t* E)
{
    cond_dist41_init(&E->cs, delta_bins * pos_bins * q3_bins * q2_bins * qual_size);
    cond_dist41_set_update_rate(&E->cs, 4);
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
    unsigned char qprev[4] = {0, 0, 0, 0};
    int delta = 6; 

    char* qs = x->qual.s;
    size_t n = x->qual.n;

    size_t i;

    for (i = 0; i < n; ++i) {
        cond_dist41_encode(E->ac, &E->cs,
                cs_index(n, i, delta,
                         charmax2(qprev[3], qprev[2]),
                         qprev[1], qprev[0]), qs[i]);


        qprev[3] = qprev[2];
        qprev[2] = qprev[1];
        qprev[1] = qprev[0];
        qprev[0] = qs[i];

        if (qprev[1] > qprev[0]) {
            delta = intmin2(delta_max - 1, delta + qprev[1] - qprev[0]);
        }
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

    unsigned char qprev[4] = {0, 0, 0, 0};
    int delta = 6; 

    size_t i;
    for (i = 0; i < n; ++i) {
        qs[i] = cond_dist41_decode(E->ac, &E->cs,
                    cs_index(n, i, delta,
                             charmax2(qprev[3], qprev[2]),
                             qprev[1], qprev[0]));

        qprev[3] = qprev[2];
        qprev[2] = qprev[1];
        qprev[1] = qprev[0];
        qprev[0] = qs[i];

        
        if (qprev[1] > qprev[0]) {
            delta = intmin2(delta_max - 1, delta + qprev[1] - qprev[0]);
        }
    }

    qual->n = n;
    qual->s[qual->n] = '\0';
}


void qualenc_reset_decoder(qualenc_t* E)
{
    ac_reset_decoder(E->ac);
}


