
#include "seqenc.h"
#include "misc.h"
#include "ac.h"
#include "dist.h"
#include <assert.h>


struct seqenc_t_
{
    /* coder */
    ac_t* ac;

    /* order of markev chain model of nucleotide probabilities */
    size_t k;

    /* bitmask for two bit encoded context */
    uint32_t ctx_mask;

    /* nucelotide probability given the last k nucleotides */
    cond_dist125_t cs;

    /* binary distribution of unique (0) / match (1) */
    dist2_t ms;

    /* distribution over match strand */
    dist2_t ss;

    /* distribution over contig number */
    // TODO

    /* distribution of edit operations, given the previous operation */
    cond_dist4_t es;
};


static void seqenc_init(seqenc_t* E, size_t k, bool decoder)
{
    E->k = k;
    size_t N = 1;

    size_t i;
    E->ctx_mask = 0;
    for (i = 0; i < k; ++i) {
        E->ctx_mask = (E->ctx_mask << 2) | 0x3;
        N *= 4;
    }

    cond_dist125_init(&E->cs, N, decoder);
    dist2_init(&E->ms, decoder);
    dist2_init(&E->ss, decoder);
    cond_dist4_init(&E->es, 16, decoder);
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
    cond_dist125_free(&E->cs);
    dist2_free(&E->ms);
    dist2_free(&E->ss);
    cond_dist4_free(&E->es);
    free(E);
}


void seqenc_encode_twobit_seq(seqenc_t* E, const twobit_t* x)
{
    dist2_encode(E->ac, &E->ms, 0);

    kmer_t u, v, w;
    size_t n = twobit_len(x);
    size_t i;
    uint32_t ctx = 0;

    /* We encode trinucleotides to improve efficiency. If the length of a read
     * is not divisible by three, we in effect pad the beginning of the read
     * with one or two Ns.
     */
    switch (n % 3) {
        case 0:
            i = 0;
            break;

        case 1:
            u = twobit_get(x, 0);
            cond_dist125_encode(E->ac, &E->cs, 0, u * 25);
            ctx = (u & 0x3) << 4;
            i = 1;
            break;

        case 2:
            v = twobit_get(x, 0);
            u = twobit_get(x, 1);
            cond_dist125_encode(E->ac, &E->cs, 0, (u * 25) + (v * 5));
            ctx = ((u & 0x3) << 4) | ((v & 0x3) << 2);
            i = 2;
            break;
    }


    for (; i < n - 2; i += 3) {
        w = twobit_get(x, i);
        v = twobit_get(x, i + 1);
        u = twobit_get(x, i + 2);
        cond_dist125_encode(E->ac, &E->cs, ctx, (u * 25) + (v * 5) + w);
        ctx = ((ctx << 6) | ((u & 0x3) << 4) | ((v & 0x3) << 2) | (w & 0x3)) & E->ctx_mask;
    }
}



void seqenc_encode_char_seq(seqenc_t* E, const char* x, size_t n)
{
    dist2_encode(E->ac, &E->ms, 0);

    kmer_t u, v, w;
    size_t i;
    uint32_t ctx = 0;

    switch (n % 3) {
        case 0:
            i = 0;
            break;

        case 1:
            u = nuc_map[(uint8_t) *x++];
            cond_dist125_encode(E->ac, &E->cs, 0, u * 25);
            ctx = (u & 0x3) << 4;
            i = 1;
            break;

        case 2:
            v = nuc_map[(uint8_t) *x++];
            u = nuc_map[(uint8_t) *x++];
            cond_dist125_encode(E->ac, &E->cs, 0, (u * 25) + (v * 5));
            ctx = ((u & 0x3) << 4) | ((v & 0x3) << 2);
            i = 2;
            break;
    }

    for (; i < n - 2; i += 3) {
        w = nuc_map[(uint8_t) *x++];
        v = nuc_map[(uint8_t) *x++];
        u = nuc_map[(uint8_t) *x++];
        cond_dist125_encode(E->ac, &E->cs, ctx, (u * 25) + (v * 5) + w);
        ctx = ((ctx << 6) | ((u & 0x3) << 4) | ((v & 0x3) << 2) | (w & 0x3)) & E->ctx_mask;
    }
}



void seqenc_encode_alignment(seqenc_t* E,
        size_t contig_idx, uint8_t strand,
        const sw_alignment_t* aln, const twobit_t* query)
{
    dist2_encode(E->ac, &E->ms, 1);
    dist2_encode(E->ac, &E->ss, strand);

    // TODO: output contig number
    // TODO: output contig offset

    kmer_t u = twobit_get(query, 0);
    uint32_t last_op = EDIT_MATCH;
    uint32_t ctx = 0;
    size_t i, j;
    for (i = 0, j = 0; i < aln->len; ++i) {
        switch (aln->ops[i]) {
            case EDIT_MATCH:
                cond_dist4_encode(E->ac, &E->es, last_op, EDIT_MATCH);

                ++j;
                u = u == 0 ? 1 : u - 1;
                ctx = (ctx << 2 | u) & E->ctx_mask;
                u = twobit_get(query, j) + 1;
                break;

            case EDIT_MISMATCH:
                cond_dist4_encode(E->ac, &E->es, last_op, EDIT_MISMATCH);
                /*cond_dist5_encode(E->ac, &E->cs, ctx, u);*/ // TODO

                ++j;
                ctx = (ctx << 2 | (u & 0x3)) & E->ctx_mask;
                u = twobit_get(query, j);
                break;

            case EDIT_Q_GAP:
                cond_dist4_encode(E->ac, &E->es, last_op, EDIT_Q_GAP);
                break;

            case EDIT_S_GAP:
                cond_dist4_encode(E->ac, &E->es, last_op, EDIT_S_GAP);
                /*cond_dist5_encode(E->ac, &E->cs, ctx, u);*/ // TODO

                ++j;
                ctx = (ctx << 2 | (u & 0x3)) & E->ctx_mask;
                u = twobit_get(query, j);
                break;
        }

        last_op = ((last_op * 4) + aln->ops[i]) % 16;
    }
}


static void seqenc_decode_seq(seqenc_t* E, seq_t* x, size_t n)
{
    kmer_t uvw, u, v, w;
    size_t i;
    uint32_t ctx = 0;

    switch (n % 3) {
        case 0:
            i = 0;
            break;

        case 1:
            u = cond_dist125_decode(E->ac, &E->cs, 0) / 25;
            x->seq.s[0] = rev_nuc_map[u];
            ctx = (u & 0x3) << 4;
            i = 1;
            break;

        case 2:
            uvw = cond_dist125_decode(E->ac, &E->cs, 0);
            u = uvw / 25;
            uvw %= 25;
            v = uvw / 5;
            x->seq.s[0] = rev_nuc_map[v];
            x->seq.s[1] = rev_nuc_map[u];
            ctx = ((u & 0x3) << 4) | ((v & 0x3) << 2);
            i = 2;
            break;
    }

    for (; i < n - 2;) {
        uvw = cond_dist125_decode(E->ac, &E->cs, ctx);
        u = uvw / 25;
        uvw %= 25;
        v = uvw / 5;
        w = uvw % 5;

        x->seq.s[i++] = rev_nuc_map[w];
        x->seq.s[i++] = rev_nuc_map[v];
        x->seq.s[i++] = rev_nuc_map[u];
        ctx = ((ctx << 6) | ((u & 0x3) << 4) | ((v & 0x3) << 2) | (w & 0x3)) & E->ctx_mask;
    }
}


void seqenc_decode(seqenc_t* E, seq_t* x, size_t n)
{
    while (n > x->seq.size) fastq_expand_str(&x->seq);

    symb_t t = dist2_decode(E->ac, &E->ms);

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


