
#include "seqenc.h"
#include "misc.h"
#include "ac.h"
#include "dist.h"
#include <assert.h>


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

static const uint8_t rev_nuc_map[5] = { 'N', 'A', 'C', 'G', 'T' };


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
    // TODO

    /* distribution of edit operations, given the previous operation */
    dist_t** es;
};


static void seqenc_init(seqenc_t* E, size_t k, bool decoder)
{
    E->k = k;
    E->N = 1;

    size_t i;
    for (i = 0; i < k; ++i) E->N *= 5;

    E->cs = malloc_or_die(E->N * sizeof(dist_t*));

    for (i = 0; i < E->N; ++i) {
        E->cs[i] = dist_alloc(5, decoder);
    }

    E->ms = dist_alloc(2, decoder);
    E->ss = dist_alloc(2, decoder);

    E->es = malloc_or_die(16 * sizeof(dist_t*));
    for (i = 0; i < 16; ++i) E->es[i] = dist_alloc(4, decoder);
}


seqenc_t* seqenc_alloc_encoder(size_t k, quip_writer_t writer, void* writer_data)
{
    seqenc_t* E = malloc_or_die(sizeof(seqenc_t));

    E->ac = ac_alloc_encoder(writer, writer_data);

    seqenc_init(E, k, false);

    return E;
}


seqenc_t* seqenc_alloc_decoder(size_t k, quip_reader_t reader, void* reader_data)
{
    seqenc_t* E = malloc_or_die(sizeof(seqenc_t));

    E->ac = ac_alloc_decoder(reader, reader_data);

    seqenc_init(E, k, true);

    return E;
}


void seqenc_free(seqenc_t* E)
{
    if (E == NULL) return;

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


static void seqenc_decode_seq(seqenc_t* E, seq_t* x, size_t n)
{
    kmer_t u;
    size_t ctx = 0;

    size_t i;
    for (i = 0; i < n; ++i) {
        u = ac_decode(E->ac, E->cs[ctx]);
        x->seq.s[i] = rev_nuc_map[u];
        ctx = (5 * ctx + u) % E->N;
    }
}


void seqenc_decode(seqenc_t* E, seq_t* x, size_t n)
{
    while (n > x->seq.size) fastq_expand_str(&x->seq);

    symb_t t = ac_decode(E->ac, E->ms);

    if (t == 0) seqenc_decode_seq(E, x, n);

    /* TODO: decoding alignments */
    assert(t == 0);
}


void seqenc_flush(seqenc_t* E)
{
    ac_flush_encoder(E->ac);
}


void seqenc_reset_decoder(seqenc_t* E)
{
    ac_reset_decoder(E->ac);
}


