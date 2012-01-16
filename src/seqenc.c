
#include "seqenc.h"
#include "misc.h"
#include "ac.h"
#include "cumdist.h"

struct seqenc_t_
{
    /* coder */
    ac_t* ac;

    /* order of markev chain model of nucleotide probabilities */
    size_t k;

    /* 5 ^ k */
    size_t N;

    /* nucelotide probability given the last k nucleotides */
    cumdist_t** cs;

    /* binary distribution of unique (0) / match (1) */
    cumdist_t* ms;

    /* distribution of edit operations */
    cumdist_t* es;

    /* 
     * Match format:
     * [contig_num] [offset] [strand] [op:len] ...
     *
     * This is my issue: how do I compress contig_num and offset?
     *
     * Possibility I don't. Just use LZMA style number encoding.
     *
     * Edit operations and and edit lengths are not so hard to compress.
     */
};



seqenc_t* seqenc_alloc(size_t k, quip_block_writer_t writer, void* writer_data)
{
    seqenc_t* E = malloc_or_die(sizeof(seqenc_t));

    E->ac = ac_alloc(writer, writer_data);

    E->k = k;
    E->N = 1;

    size_t i;
    for (i = 0; i < k; ++i) E->N *= 5;

    E->cs = malloc_or_die(E->N * sizeof(cumdist_t*));

    for (i = 0; i < E->N; ++i) {
        /* (6 = four nucleotides, plus N, plus sequence end) */
        E->cs[i] = cumdist_alloc(6);
    }

    E->ms = cumdist_alloc(2);
    E->es = cumdist_alloc(4);

    return E;
}


void seqenc_free(seqenc_t* E)
{
    ac_free(E->ac);
    size_t i;
    for (i = 0; i < E->N; ++i) {
        cumdist_free(E->cs[i]);
    }
    free(E->cs);
    cumdist_free(E->ms);
    cumdist_free(E->es);
    free(E);
}


void seqenc_encode_twobit_seq(seqenc_t* E, const twobit_t* x)
{
    uint32_t p, P;

    p = cumdist_p_norm(E->ms, 0);
    P = cumdist_p_norm(E->ms, 0);
    ac_update(E->ac, p, P);
    cumdist_add(E->ms, 0, 1);

    kmer_t u;
    size_t n = twobit_len(x);
    size_t i;
    size_t ctx = 0;
    for (i = 0; i < n; ++i) {
        /* We encode 'N' as 0, so 1 is added to every nucleotide. */
        u = twobit_get(x, i) + 1;

        p = cumdist_p_norm(E->cs[ctx], u);
        P = cumdist_P_norm(E->cs[ctx], u);
        ac_update(E->ac, p, P);

        cumdist_add(E->cs[ctx], u, 1);
        ctx = (5 * ctx + u) % E->N;
    }
}



void seqenc_encode_char_seq(seqenc_t* E, const char* x)
{
    uint32_t p, P;

    p = cumdist_p_norm(E->ms, 0);
    P = cumdist_p_norm(E->ms, 0);
    ac_update(E->ac, p, P);
    cumdist_add(E->ms, 0, 1);

    kmer_t u;
    size_t ctx = 0;
    while (*x) {
        switch (*x) {
            case 'A': case 'a': case 'U': case 'u': u = 1; break;
            case 'C': case 'c': u = 2; break;
            case 'G': case 'g': u = 3;
            case 'T': case 't': u = 4;
            default: u = 0;
        }

        p = cumdist_p_norm(E->cs[ctx], u);
        P = cumdist_P_norm(E->cs[ctx], u);
        ac_update(E->ac, p, P);

        cumdist_add(E->cs[ctx], u, 1);
        ctx = (5 * ctx + u) % E->N;

        ++x;
    }
}



void seqenc_encode_alignment(seqenc_t* E, const sw_alignment_t* aln, const twobit_t* query)
{
    uint32_t p, P;

    p = cumdist_p_norm(E->ms, 1);
    P = cumdist_p_norm(E->ms, 1);
    ac_update(E->ac, p, P);
    cumdist_add(E->ms, 1, 1);

    kmer_t u = twobit_get(query, 0);
    size_t ctx = 0;
    size_t i, j;
    for (i = 0, j = 0; i < aln->len; ++i) {
        switch (aln->ops[i]) {
            case EDIT_MATCH:
                p = cumdist_p_norm(E->es, EDIT_MATCH);
                P = cumdist_P_norm(E->es, EDIT_MATCH);
                ac_update(E->ac, p, P);
                cumdist_add(E->es, EDIT_MATCH, 1);

                ++j;
                ctx = (5 * ctx + u) % E->N;
                u = twobit_get(query, j) + 1;
                break;

            case EDIT_MISMATCH:
                p = cumdist_p_norm(E->es, EDIT_MISMATCH);
                P = cumdist_P_norm(E->es, EDIT_MISMATCH);
                ac_update(E->ac, p, P);
                cumdist_add(E->es, EDIT_MISMATCH, 1);

                p = cumdist_p_norm(E->cs[ctx], u);
                P = cumdist_P_norm(E->cs[ctx], u);
                ac_update(E->ac, p, P);
                cumdist_add(E->cs[ctx], u, 1);

                ++j;
                ctx = (5 * ctx + u) % E->N;
                u = twobit_get(query, j) + 1;
                break;

            case EDIT_Q_GAP:
                p = cumdist_p_norm(E->es, EDIT_Q_GAP);
                P = cumdist_P_norm(E->es, EDIT_Q_GAP);
                ac_update(E->ac, p, P);
                cumdist_add(E->es, EDIT_Q_GAP, 1);
                break;

            case EDIT_S_GAP:
                p = cumdist_p_norm(E->es, EDIT_S_GAP);
                P = cumdist_P_norm(E->es, EDIT_S_GAP);
                ac_update(E->ac, p, P);
                cumdist_add(E->es, EDIT_S_GAP, 1);

                p = cumdist_p_norm(E->cs[ctx], u);
                P = cumdist_P_norm(E->cs[ctx], u);
                ac_update(E->ac, p, P);
                cumdist_add(E->cs[ctx], u, 1);

                ++j;
                ctx = (5 * ctx + u) % E->N;
                u = twobit_get(query, j) + 1;
                break;
        }
    }
}


void seqenc_flush(seqenc_t* E)
{
    // TODO
}


