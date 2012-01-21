
#include "seqenc.h"
#include "misc.h"
#include "ac.h"
#include "dist.h"


/* map nucleotide ascii characters to numbers */
static const uint8_t nuc_map[256] =
  { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 1, 0, 2, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 1, 0, 2, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

struct seqenc_t_
{
    /* coder */
    ac_t* ac;

    /* order of markev chain model of nucleotide probabilities */
    size_t k;

    /* 5 ^ k */
    size_t N;

    /* nucelotide probability given the last k nucleotides */
    dist_t** cs;

    /* binary distribution of unique (0) / match (1) */
    dist_t* ms;

    /* distribution over match strand */
    dist_t* ss;

    /* distribution over contig number */
    // XXX: fuck!!!! How do I do this? 

    /* distribution of edit operations, given the previous operation */
    dist_t** es;
};



seqenc_t* seqenc_alloc(size_t k, quip_writer_t writer, void* writer_data)
{
    seqenc_t* E = malloc_or_die(sizeof(seqenc_t));

    E->ac = ac_alloc_encoder(writer, writer_data);

    E->k = k;
    E->N = 1;

    size_t i;
    for (i = 0; i < k; ++i) E->N *= 5;

    E->cs = malloc_or_die(E->N * sizeof(dist_t*));

    for (i = 0; i < E->N; ++i) {
        E->cs[i] = dist_alloc_encode(5);
    }

    E->ms = dist_alloc_encode(2);
    E->ss = dist_alloc_encode(2);

    E->es = malloc_or_die(16 * sizeof(dist_t*));
    for (i = 0; i < 16; ++i) E->es[i] = dist_alloc_encode(4);

    return E;
}


void seqenc_free(seqenc_t* E)
{
    ac_free(E->ac);
    size_t i;
    for (i = 0; i < E->N; ++i) {
        dist_free(E->cs[i]);
    }
    free(E->cs);
    dist_free(E->ms);
    dist_free(E->ss);
    for (i = 0; i < 16; ++i) {
        dist_free(E->es[i]);
    }
    free(E->es);
    free(E);
}


void seqenc_encode_twobit_seq(seqenc_t* E, const twobit_t* x)
{
    ac_encode(E->ac, E->ms, 0);

    kmer_t u;
    size_t n = twobit_len(x);
    size_t i;
    size_t ctx = 0;
    for (i = 0; i < n; ++i) {
        /* We encode 'N' as 0, so 1 is added to every nucleotide. */
        u = twobit_get(x, i) + 1;
        ac_encode(E->ac, E->cs[ctx], u);
        ctx = (5 * ctx + u) % E->N;
    }
}



void seqenc_encode_char_seq(seqenc_t* E, const char* x)
{
    ac_encode(E->ac, E->ms, 0);

    kmer_t u;
    size_t ctx = 0;
    while (*x) {
        u = nuc_map[(uint8_t) *x];
        ac_encode(E->ac, E->cs[ctx], u);
        ctx = (5 * ctx + u) % E->N;
        ++x;
    }
}



void seqenc_encode_alignment(seqenc_t* E,
        size_t contig_idx, uint8_t strand,
        const sw_alignment_t* aln, const twobit_t* query)
{
    ac_encode(E->ac, E->ms, 1);
    ac_encode(E->ac, E->ss, strand);

    // TODO: output contig number
    // TODO: output contig offset

    kmer_t u = twobit_get(query, 0);
    size_t last_op = EDIT_MATCH;
    size_t ctx = 0;
    size_t i, j;
    for (i = 0, j = 0; i < aln->len; ++i) {
        switch (aln->ops[i]) {
            case EDIT_MATCH:
                ac_encode(E->ac, E->es[last_op], EDIT_MATCH);

                ++j;
                ctx = (5 * ctx + u) % E->N;
                u = twobit_get(query, j) + 1;
                break;

            case EDIT_MISMATCH:
                ac_encode(E->ac, E->es[last_op], EDIT_MISMATCH);
                ac_encode(E->ac, E->cs[ctx], u);

                ++j;
                ctx = (5 * ctx + u) % E->N;
                u = twobit_get(query, j) + 1;
                break;

            case EDIT_Q_GAP:
                ac_encode(E->ac, E->es[last_op], EDIT_Q_GAP);
                break;

            case EDIT_S_GAP:
                ac_encode(E->ac, E->es[last_op], EDIT_S_GAP);
                ac_encode(E->ac, E->cs[ctx], u);

                ++j;
                ctx = (5 * ctx + u) % E->N;
                u = twobit_get(query, j) + 1;
                break;
        }

        last_op = ((last_op * 4) + aln->ops[i]) % 16;
    }
}


void seqenc_flush(seqenc_t* E)
{
    ac_flush_encoder(E->ac);
}


