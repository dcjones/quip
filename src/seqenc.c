
#include "seqenc.h"
#include "misc.h"
#include "ac.h"
#include "dist.h"
#include <assert.h>
#include <string.h>


/* map nucleotide ascii characters to numbers */
static const uint8_t nuc_map[256] =
  { 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 0, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 0, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 0 };

static const uint8_t rev_nuc_map[5] = { 'A', 'C', 'G', 'T', 'N' };

static const size_t insert_nuc_k = 4;


struct seqenc_t_
{
    /* coder */
    ac_t* ac;

    /* order of markev chain model of nucleotide probabilities */
    size_t k;

    /* bitmask for two bit encoded context */
    uint32_t ctx_mask;

    /* nucelotide probability given the last k nucleotides */
    cond_dist16_t cs;

    /* binary distribution of unique (0) / match (1) */
    dist2_t ms;

    /* distribution over match strand */
    dist2_t ss;

    /* distribution over contig number */
    dist4_t        d_contig_idx_bytes;
    cond_dist256_t d_contig_idx;

    /* distribution over the offset into the contig */
    dist4_t        d_contig_off_bytes;
    cond_dist256_t d_contig_off;

    /* distribution over inserted nucleotides */
    cond_dist4_t d_ins_nuc;
    uint32_t ins_ctx_mask;

    /* distribution of edit operations, given the previous operation */
    cond_dist4_t es;

    /* store contigs, for decoding */
    twobit_t** contigs;
    size_t contig_count;
    size_t remaining_contigs;
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

    cond_dist16_init(&E->cs, N, decoder);
    dist2_init(&E->ms, decoder);
    dist2_init(&E->ss, decoder);

    dist4_init(&E->d_contig_idx_bytes, decoder);
    cond_dist256_init(&E->d_contig_idx, 4, decoder);

    dist4_init(&E->d_contig_off_bytes, decoder);
    cond_dist256_init(&E->d_contig_off, 4, decoder);

    N = 1;
    E->ins_ctx_mask = 0;
    for (i = 0; i < insert_nuc_k; ++i) {
        E->ins_ctx_mask = (E->ins_ctx_mask << 2) | 0x3;
        N *= 4;
    }
    cond_dist4_init(&E->d_ins_nuc, N, decoder);

    cond_dist4_init(&E->es, 16, decoder);

    E->contigs = NULL;
    E->contig_count = 0;
    E->remaining_contigs = 0;
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


void seqenc_prepare_decoder(seqenc_t* E, size_t contig_count)
{
    size_t i;
    for (i = 0; i < E->contig_count; ++i) {
        twobit_free(E->contigs[i]);
    }
    free(E->contigs);

    E->contig_count = 0;
    E->remaining_contigs = contig_count;
    E->contigs = malloc_or_die(E->remaining_contigs * sizeof(twobit_t*));
}


void seqenc_free(seqenc_t* E)
{
    if (E == NULL) return;

    ac_free(E->ac);
    cond_dist16_free(&E->cs);
    dist2_free(&E->ms);
    dist2_free(&E->ss);

    dist4_free(&E->d_contig_idx_bytes);
    cond_dist256_free(&E->d_contig_idx);

    dist4_free(&E->d_contig_off_bytes);
    cond_dist256_free(&E->d_contig_off);

    cond_dist4_free(&E->d_ins_nuc);

    cond_dist4_free(&E->es);

    size_t i;
    for (i = 0; i < E->contig_count; ++i) {
        twobit_free(E->contigs[i]);
    }
    free(E->contigs);


    free(E);
}


void seqenc_encode_twobit_seq(seqenc_t* E, const twobit_t* x)
{
    dist2_encode(E->ac, &E->ms, 0);

    kmer_t uv;
    size_t n = twobit_len(x);
    size_t i;
    uint32_t ctx = 0;
    for (i = 0; i < n - 1; i += 2) {
        uv = (twobit_get(x, i) << 2) | twobit_get(x, i + 1);
        cond_dist16_encode(E->ac, &E->cs, ctx, uv);
        ctx = ((ctx << 4) | uv) & E->ctx_mask;
    }

    /* handle odd read lengths */
    if (i < n) {
        uv = twobit_get(x, i);
        cond_dist16_encode(E->ac, &E->cs, ctx, uv);
    }
}



void seqenc_encode_char_seq(seqenc_t* E, const char* x, size_t len)
{
    dist2_encode(E->ac, &E->ms, 0);

    kmer_t uv;
    uint32_t ctx = 0;
    size_t i = 0;
    size_t j = 0;

    for (i = 0; i < len; ++i) {
        if (x[i] == 'N') break;
    }

    /* the sequence contains Ns, encode around them */
    if (i < len) {
        i = 0;
        while (i < len - 1) {
            while (i < len && x[i] == 'N') ++i;
            if (i == len) break;

            j = i + 1;
            while (j < len && x[j] == 'N') ++j;
            if (j == len) break;

            uv = (nuc_map[(uint8_t) x[i]] << 2) | nuc_map[(uint8_t) x[j]];
            cond_dist16_encode(E->ac, &E->cs, ctx, uv);
            ctx = ((ctx << 4) | uv) & E->ctx_mask;

            i = j + 1;
        }

        /* handle odd read lengths */
        if (i < len && x[i] != 'N') {
            uv = nuc_map[(uint8_t) x[i]];
            cond_dist16_encode(E->ac, &E->cs, ctx, uv);
        }
    }
    /* no Ns, use a slightly more efficent loop */
    else {
        for (i = 0; i < len - 1; i += 2) {
            uv = (nuc_map[(uint8_t) x[i]] << 2) | nuc_map[(uint8_t) x[i + 1]];
            cond_dist16_encode(E->ac, &E->cs, ctx, uv);
            ctx = ((ctx << 4) | uv) & E->ctx_mask;
        }

        /* handle odd read lengths */
        if (i == len - 1) {
            uv = nuc_map[(uint8_t) x[i]];
            cond_dist16_encode(E->ac, &E->cs, ctx, uv);
        }
    }
}



void seqenc_encode_alignment(seqenc_t* E,
        size_t contig_idx, uint8_t strand,
        const sw_alignment_t* aln, const twobit_t* query)
{
    dist2_encode(E->ac, &E->ms, 1);

    uint32_t bytes;
    uint32_t z;


    /* encode contig index */
    bytes = 0;
    z = contig_idx;
    while (z > 0) {
        z >>= 8;
        ++bytes;
    }
    dist4_encode(E->ac, &E->d_contig_idx_bytes, bytes);

    z = contig_idx;
    while (z > 0) {
        cond_dist256_encode(E->ac, &E->d_contig_idx, bytes - 1, z & 0xff);
        z >>= 8;
    }


    /* encode strand */
    dist2_encode(E->ac, &E->ss, strand);


    /* encode contig offset */
    bytes = 0;
    z = aln->spos;
    while (z > 0) {
        z >>= 8;
        ++bytes;
    }
    dist4_encode(E->ac, &E->d_contig_idx_bytes, bytes);

    z = aln->spos;
    while (z > 0) {
        cond_dist256_encode(E->ac, &E->d_contig_idx, bytes - 1, z & 0xff);
        z >>= 8;
    }


    /* encode edit operations */
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
                cond_dist4_encode(E->ac, &E->d_ins_nuc, ctx, u);

                ++j;
                ctx = (ctx << 2 | (u & 0x3)) & E->ctx_mask;
                u = twobit_get(query, j);
                break;

            case EDIT_Q_GAP:
                cond_dist4_encode(E->ac, &E->es, last_op, EDIT_Q_GAP);
                break;

            case EDIT_S_GAP:
                cond_dist4_encode(E->ac, &E->es, last_op, EDIT_S_GAP);
                cond_dist4_encode(E->ac, &E->d_ins_nuc, ctx, u);

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
    kmer_t uv, u, v;
    uint32_t ctx = 0;
    size_t i, j;

    const char* qs = x->qual.s;

    for (i = 0; i < n; ++i) {
        if (qs[i] == 0) break;
    }

    /* sequence contains Ns */
    if (i < n) {

        memset(x->seq.s, 'N', n * sizeof(char));

        i = 0;
        while (i < n - 1) {
            while (i < n && qs[i] == 0) ++i;
            if (i == n) break;

            j = i + 1;
            while (j < n && qs[j] == 0) ++j;
            if (j == n) break;

            uv = cond_dist16_decode(E->ac, &E->cs, ctx);
            u = uv >> 2;
            v = uv & 0x3;

            x->seq.s[i] = rev_nuc_map[u];
            x->seq.s[j] = rev_nuc_map[v];

            ctx = ((ctx << 4) | uv) & E->ctx_mask;

            i = j + 1;
        }

        if (i < n && qs[i] != 0) {
            uv = cond_dist16_decode(E->ac, &E->cs, ctx);
            v = uv & 0x3;
            x->seq.s[i] = rev_nuc_map[v];
        }
    }
    else {
        for (i = 0; i < n - 1;) {
            uv = cond_dist16_decode(E->ac, &E->cs, ctx);
            u = uv >> 2;
            v = uv & 0x3;
            x->seq.s[i++] = rev_nuc_map[u];
            x->seq.s[i++] = rev_nuc_map[v];
            ctx = ((ctx << 4) | uv) & E->ctx_mask;
        }
    }
}


static void seqenc_decode_alignment(seqenc_t* E, seq_t* x, size_t n)
{
    uint32_t contig_idx, spos;
    uint32_t bytes;
    uint32_t b;
    uint32_t z;

    /* decode contig index */
    bytes = dist4_decode(E->ac, &E->d_contig_idx_bytes);
    z = 0;
    for (b = 0; b < bytes; ++b) {
        z |= cond_dist256_decode(E->ac, &E->d_contig_idx, bytes - b - 1) << (8 * b);
    }
    contig_idx = z;


    /* decode strand */
    uint8_t strand = dist2_decode(E->ac, &E->ss);


    /* decode offset */
    bytes = dist4_decode(E->ac, &E->d_contig_off_bytes);
    z = 0;
    for (b = 0; b < bytes; ++b) {
        z |= cond_dist256_decode(E->ac, &E->d_contig_idx, bytes - b - 1) << (8 * b);
    }
    spos = z;


    /* decode edit operations */
    edit_op_t op, last_op = EDIT_MATCH;
    uint32_t ctx;

    size_t i = 0;
    while (i < n) {
        op = cond_dist4_decode(E->ac, &E->es, last_op);

        switch (op) {
            case EDIT_MATCH:
                // TODO
                break;

            case EDIT_MISMATCH:
                // TODO
                break;

            case EDIT_Q_GAP:
                // TODO
                break;

            case EDIT_S_GAP:
                // TODO
                break;
        }

        last_op = ((last_op * 4) + op) % 16;
    }
}



bool seqenc_decode(seqenc_t* E, seq_t* x, size_t n)
{
    while (n > x->seq.size) fastq_expand_str(&x->seq);

    symb_t t = dist2_decode(E->ac, &E->ms);

    if (t == 0) seqenc_decode_seq(E, x, n);
    else        seqenc_decode_alignment(E, x, n);

    if (E->remaining_contigs > 0) {
        E->contigs[E->contig_count] = twobit_alloc_n(x->seq.n);
        twobit_copy_n(E->contigs[E->contig_count], x->seq.s, x->seq.n);
        E->contig_count++;
        E->remaining_contigs--;
        return false;
    }

    return true;
}


void seqenc_flush(seqenc_t* E)
{
    ac_flush_encoder(E->ac);
}


void seqenc_reset_decoder(seqenc_t* E)
{
    size_t i;
    for (i = 0; i < E->contig_count; ++i) {
        twobit_free(E->contigs[i]);
    }
    E->contig_count = 0;
    E->remaining_contigs = 0;
    free(E->contigs);
    E->contigs = NULL;

    ac_reset_decoder(E->ac);
}


