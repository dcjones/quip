
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

    /* binary distribution of unique / match */
    cumdist_t* ms;

    /* distribution of edit operations */
    cumdist_t* es;

    /* distribution over edit operation lengths */
    cumdist_t* ls;

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
    /*E->es = cumdist_alloc();*/

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
    free(E);
}


void seqenc_encode_twobit_seq(seqenc_t* E, const twobit_t* x)
{
    kmer_t u;
    size_t n = twobit_len(x);
    uint32_t p, P, Z;
    size_t i;
    size_t ctx = 0;
    for (i = 0; i < n; ++i) {
        /* We encode 'N' as 0, so 1 is added to every nucleotide. */
        u = twobit_get(x, i) + 1;

        Z = cumdist_Z(E->cs[ctx]);
        p = ((uint64_t) cumdist_p(E->cs[ctx], u) << 32) / Z;
        P = ((uint64_t) cumdist_P(E->cs[ctx], u) << 32) / Z;
        if (P == 0) P = 0xffffffff;
        ac_update(E->ac, p, P);

        cumdist_add(E->cs[ctx], u, 1);
        ctx = (5 * ctx + u) % E->N;
    }

    /* encode end-of-sequence sigil */
    Z = cumdist_Z(E->cs[ctx]);
    p = ((uint64_t) cumdist_p(E->cs[ctx], 5) << 32) / Z;
    P = ((uint64_t) cumdist_P(E->cs[ctx], 5) << 32) / Z;
    if (P == 0) P = 0xffffffff;
    ac_update(E->ac, p, P);

    cumdist_add(E->cs[ctx], 5, 1);
}



void seqenc_encode_char_seq(seqenc_t* E, const char* x)
{
    kmer_t u;
    uint32_t p, P;
    size_t ctx = 0;
    while (*x) {
        switch (*x) {
            case 'A': case 'a': case 'U': case 'u': u = 1; break;
            case 'C': case 'c': u = 2; break;
            case 'G': case 'g': u = 3;
            case 'T': case 't': u = 4;
            default: u = 0;
        }

        p = cumdist_p(E->cs[ctx], u);
        P = cumdist_P(E->cs[ctx], u);
        ac_update(E->ac, p, P);

        cumdist_add(E->cs[ctx], u, 1);
        ctx = (5 * ctx + u) % E->N;

        ++x;
    }

    /* encode end-of-sequence sigil */
    p = cumdist_p(E->cs[ctx], 5);
    P = cumdist_P(E->cs[ctx], 5);
    ac_update(E->ac, p, P);

    cumdist_add(E->cs[ctx], 5, 1);
}



