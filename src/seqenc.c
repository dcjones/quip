
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
    cumdist_t** C;

    /* distribution of edit operations */
    cumdist_t* D;

    /* distribution over edit operation lengths */
    cumdist_t* L;



    /* This follows lzma somewhat.
     * The poosible operations are essetionally,
     * LZ-style dictionary match, or encoded nucleotide.
     *
     *
     * NOMATCH : 
     * MATCH   : CONTIG_NUM : OFFSET : STRAND : LEN
     *
     * But this is not quite right. We only allow one alignment, so there is no
     * point in repeatedly specifying the contig number offset and strand.
     *
     * We do not necessarily need to compress the
     * contig_num : offset : strand
     *
     *
     * So the whole format will look like this:
     *
     *  [contig count][contigs][reads]
     *
     *  Each sequence is encoded as follows:
     *
     *  Matches:
     *    [match][contig_num][offset][strand]
     *    followed by a sequence of edit operations (i.e. match, nonmatch)
     *
     *
     * Maybe reads should be prefixed by their lengths. We have to include in
     * the header the distribution of lengths, but this should be relatively
     * tight.
     *
     *
     * What about collapsing identical reads??
     *
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

    E->C = malloc_or_die(sizeof(cumdist_t*));

    for (i = 0; i < E->N; ++i) {
        /* (6 = four nucleotides, plus N, plus sequence end) */
        E->C[i] = cumdist_alloc(6);
    }

    return E;
}


void seqenc_free(seqenc_t* E)
{
    ac_free(E->ac);
    size_t i;
    for (i = 0; i < E->N; ++i) {
        cumdist_free(E->C[i]);
    }
    free(E->C);
    free(E);
}


void seqenc_encode_twobit_seq(seqenc_t* E, const twobit_t* x)
{
    kmer_t u;
    size_t n = twobit_len(x);
    uint32_t p, P;
    size_t i;
    size_t ctx = 0;
    for (i = 0; i < n; ++i) {
        /* We encode 'N' as 0, so 1 is added to every nucleotide. */
        u = twobit_get(x, i) + 1;

        p = cumdist_p(E->C[ctx], u);
        P = cumdist_P(E->C[ctx], u);
        ac_update(E->ac, p, P);

        cumdist_add(E->C[ctx], u, 1);
        ctx = (5 * ctx + u) % E->N;
    }

    /* encode end-of-sequence sigil */
    p = cumdist_p(E->C[ctx], 5);
    P = cumdist_P(E->C[ctx], 5);
    ac_update(E->ac, p, P);

    cumdist_add(E->C[ctx], 5, 1);
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

        p = cumdist_p(E->C[ctx], u);
        P = cumdist_P(E->C[ctx], u);
        ac_update(E->ac, p, P);

        cumdist_add(E->C[ctx], u, 1);
        ctx = (5 * ctx + u) % E->N;

        ++x;
    }

    /* encode end-of-sequence sigil */
    p = cumdist_p(E->C[ctx], 5);
    P = cumdist_P(E->C[ctx], 5);
    ac_update(E->ac, p, P);

    cumdist_add(E->C[ctx], 5, 1);
}



