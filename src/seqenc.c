
#include "seqenc.h"
#include "seqenc_prior.h"
#include "misc.h"
#include "ac.h"
#include "dist.h"
#include <assert.h>
#include <string.h>


/* Number of mismatchisg reads before a contig position is flipped. */
const uint16_t mismatch_patch_factor = 5;
const uint16_t mismatch_patch_cutoff = 5;

/* Order of the markov chain assigning probabilities to dinucleotides. */
static const size_t k = 11;

/* Use a seperate model for the first n dinucleotides. This is primarily to
 * account for positional sequence bias that is sommon in short read sequencing.  */
static const size_t prefix_len = 4;

/* Order of the markov chain assigning probabilities to inserted nucleotides in
 * alignments. */
static const size_t insert_nuc_k = 4;

/* order of the markov-chain assigning probabilities to edit operations in
 * encoded alignment.s
 */
static const size_t edit_op_k = 6;

/* The rate at which the nucleotide markov chain is updated. */
static const size_t seq_update_rate = 2;



struct seqenc_t_
{
    /* coder */
    ac_t* ac;

    /* bitmask for two bit encoded context */
    uint32_t ctx_mask;
    uint32_t ctx0_mask;

    /* nucelotide probability given the last k nucleotides */
    cond_dist16_t cs;

    /* special case models for the first 2 * prefix_les positions. */
    cond_dist16_t cs0[5];

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

    /* distribution of edit operations, given the previous operation and
     * position within the read */
    cond_dist4_t d_edit_op;
    uint32_t edit_op_mask;

    /* contig set, used in calls to seqenc_decode_alignment */
    twobit_t** contigs;

    /* number of contigs read so far */
    uint32_t contig_count;

    /* how many of the leading sequences are contigs */
    uint32_t expected_contigs;

    /* lengths of the expected contigs */
    uint32_t* contig_lens;

    /* tallies of mismatches */
    uint16_t** mismatch_tally;
};


static void seqenc_setprior(seqenc_t* E)
{
    size_t N = 1 << (2 * seqenc_prior_k);
    size_t M = 1 << (2 * (k - seqenc_prior_k));
    size_t u;
    size_t i;

    for (u = 0; u < N; ++u) {
        for (i = 0; i < 16; ++i) {
            E->cs.xss[u].xs[i].count = seqenc_prior[(u << 4) + i];
        }
        dist16_update(&E->cs.xss[u]);

        for (i = 1; i < M; ++i) {
            memcpy(&E->cs.xss[(i << (2 * seqenc_prior_k)) | u],
                   &E->cs.xss[u], sizeof(dist16_t));
        }
    }
}



static void seqenc_init(seqenc_t* E)
{
    size_t N = 1 << (2 * k);
    E->ctx_mask = N - 1;

    cond_dist16_init(&E->cs, N);
    cond_dist16_set_update_rate(&E->cs, seq_update_rate);

    size_t i;
    for (i = 0; i < prefix_len; ++i) {
        cond_dist16_init(&E->cs0[i], 1 << (4 * i));
        cond_dist16_set_update_rate(&E->cs0[i], seq_update_rate);
    }

    /* choose an initial distribution that is slightly more informed than
     * uniform */
    seqenc_setprior(E);

    dist2_init(&E->ms);
    dist2_init(&E->ss);

    dist4_init(&E->d_contig_idx_bytes);
    cond_dist256_init(&E->d_contig_idx, 4);

    dist4_init(&E->d_contig_off_bytes);
    cond_dist256_init(&E->d_contig_off, 4);


    N = 1 << (2 * insert_nuc_k);
    E->ins_ctx_mask = N -1;

    cond_dist4_init(&E->d_ins_nuc, N);


    N = 1 << (2 * edit_op_k);
    E->edit_op_mask = N - 1;

    cond_dist4_init(&E->d_edit_op, N);

    E->contigs          = NULL;
    E->contig_lens      = NULL;
    E->contig_count     = 0;
    E->expected_contigs = 0;
    E->mismatch_tally   = NULL;
}


seqenc_t* seqenc_alloc_encoder(quip_writer_t writer, void* writer_data)
{
    seqenc_t* E = malloc_or_die(sizeof(seqenc_t));

    E->ac = ac_alloc_encoder(writer, writer_data);

    seqenc_init(E);

    return E;
}


seqenc_t* seqenc_alloc_decoder(quip_reader_t reader, void* reader_data)
{
    seqenc_t* E = malloc_or_die(sizeof(seqenc_t));

    E->ac = ac_alloc_decoder(reader, reader_data);

    seqenc_init(E);

    return E;
}


void seqenc_free(seqenc_t* E)
{
    if (E == NULL) return;

    ac_free(E->ac);
    cond_dist16_free(&E->cs);

    size_t i;
    for (i = 0; i < prefix_len; ++i) {
        cond_dist16_free(&E->cs0[i]);
    }

    cond_dist256_free(&E->d_contig_idx);
    cond_dist256_free(&E->d_contig_off);

    cond_dist4_free(&E->d_ins_nuc);
    cond_dist4_free(&E->d_edit_op);

    for (i = 0; i < E->contig_count; ++i) {
        twobit_free(E->contigs[i]);
        free(E->mismatch_tally[i]);
    }
    free(E->contigs);
    free(E->contig_lens);
    free(E->mismatch_tally);

    free(E);
}



void seqenc_encode_twobit_seq(seqenc_t* E, const twobit_t* x)
{
    dist2_encode(E->ac, &E->ms, 0);

    kmer_t uv;
    size_t n = twobit_len(x);
    size_t i;
    uint32_t ctx = 0;

    for (i = 0; i < n - 1 && i / 2 < prefix_len; i += 2) {
        uv = (twobit_get(x, i) << 2) | twobit_get(x, i + 1);
        cond_dist16_encode(E->ac, &E->cs0[i/2], ctx, uv);
        ctx = ((ctx << 4) | uv) & E->ctx_mask;
    }

    for (; i < n - 1; i += 2) {
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

            uv = (chartokmer[(uint8_t) x[i]] << 2) | chartokmer[(uint8_t) x[j]];
            cond_dist16_encode(E->ac, &E->cs, ctx, uv);
            ctx = ((ctx << 4) | uv) & E->ctx_mask;

            i = j + 1;
        }

        /* handle odd read lengths */
        if (i < len && x[i] != 'N') {
            uv = chartokmer[(uint8_t) x[i]];
            cond_dist16_encode(E->ac, &E->cs, ctx, uv);
        }
    }
    /* no Ns, use a slightly more efficent loop */
    else {
        for (i = 0; i < len - 1 && i / 2 < prefix_len; i += 2) {
            uv = (chartokmer[(uint8_t) x[i]] << 2) | chartokmer[(uint8_t) x[i + 1]];
            cond_dist16_encode(E->ac, &E->cs0[i/2], ctx, uv);
            ctx = ((ctx << 4) | uv) & E->ctx_mask;
        }

        for (; i < len - 1; i += 2) {
            uv = (chartokmer[(uint8_t) x[i]] << 2) | chartokmer[(uint8_t) x[i + 1]];
            cond_dist16_encode(E->ac, &E->cs, ctx, uv);
            ctx = ((ctx << 4) | uv) & E->ctx_mask;
        }

        /* handle odd read lengths */
        if (i == len - 1) {
            uv = chartokmer[(uint8_t) x[i]];
            cond_dist16_encode(E->ac, &E->cs, ctx, uv);
        }
    }
}



void seqenc_encode_alignment(seqenc_t* E,
        size_t contig_idx, uint8_t strand,
        const sw_alignment_t* aln, const twobit_t* query)
{
    /*print_alignment(stderr, aln);*/

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
    if (bytes == 0) bytes = 1;

    dist4_encode(E->ac, &E->d_contig_idx_bytes, bytes - 1);

    z = contig_idx;
    while (bytes--) {
        cond_dist256_encode(E->ac, &E->d_contig_idx, bytes, z & 0xff);
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
    if (bytes == 0) bytes = 1;

    dist4_encode(E->ac, &E->d_contig_off_bytes, bytes - 1);

    z = aln->spos;
    while (bytes--) {
        cond_dist256_encode(E->ac, &E->d_contig_off, bytes, z & 0xff);
        z >>= 8;
    }


    /* encode edit operations */
    kmer_t u = twobit_get(query, 0);
    uint32_t last_op = EDIT_MATCH;
    kmer_t ctx = 0;
    size_t i; /* position within the alignment */
    size_t j; /* position within the read sequence */
    for (i = 0, j = 0; i < aln->len; ++i) {
        switch (aln->ops[i]) {
            case EDIT_MATCH:
                cond_dist4_encode(E->ac, &E->d_edit_op, last_op, EDIT_MATCH);

                ++j;
                ctx = ((ctx << 2) | u) & E->ins_ctx_mask;
                u = twobit_get(query, j);
                break;

            case EDIT_MISMATCH:
                cond_dist4_encode(E->ac, &E->d_edit_op, last_op, EDIT_MISMATCH);
                cond_dist4_encode(E->ac, &E->d_ins_nuc, ctx, u);

                ++j;
                ctx = ((ctx << 2) | u) & E->ins_ctx_mask;
                u = twobit_get(query, j);
                break;

            case EDIT_Q_GAP:
                cond_dist4_encode(E->ac, &E->d_edit_op, last_op, EDIT_Q_GAP);
                break;

            case EDIT_S_GAP:
                cond_dist4_encode(E->ac, &E->d_edit_op, last_op, EDIT_S_GAP);
                cond_dist4_encode(E->ac, &E->d_ins_nuc, ctx, u);

                ++j;
                ctx = ((ctx << 2) | u) & E->ins_ctx_mask;
                u = twobit_get(query, j);
                break;
        }

        last_op = ((last_op << 2) | aln->ops[i]) & E->edit_op_mask;
    }

    assert(j == twobit_len(query));
}


void seqenc_prepare_decoder(seqenc_t* E, uint32_t n, const uint32_t* lens)
{
    size_t i;
    for (i = 0; i < E->contig_count; ++i) {
        twobit_free(E->contigs[i]);
    }
    free(E->contigs);
    E->contigs = NULL;

    free(E->contig_lens);
    E->contig_lens = NULL;

    E->contig_count = 0;
    E->expected_contigs = n;

    E->contigs     = malloc_or_die(n * sizeof(twobit_t*));
    memset(E->contigs, 0, n * sizeof(twobit_t*));

    E->contig_lens = malloc_or_die(n * sizeof(uint32_t));

    for (i = 0; i < n; ++i) E->contig_lens[i] = lens[i];

    E->mismatch_tally = malloc_or_die(n * sizeof(uint16_t*));
    for (i = 0; i < n; ++i) {
        E->mismatch_tally[i] = malloc_or_die(4 * E->contig_lens[i] * sizeof(uint16_t));
        memset(E->mismatch_tally[i], 0, 4 * E->contig_lens[i] * sizeof(uint16_t));
    }
}



static void seqenc_decode_seq(seqenc_t* E, seq_t* x, size_t n)
{
    while (n >= x->seq.size) fastq_expand_str(&x->seq);

    kmer_t uv, u, v;
    uint32_t ctx = 0;
    size_t i, j;

    const char* qs = x->qual.s;

    for (i = 0; i < x->qual.n; ++i) {
        if (qs[i] == 0) break;
    }

    /* sequence contains Ns */
    if (i < n && x->qual.n > 0) {

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

            x->seq.s[i] = kmertochar[u];
            x->seq.s[j] = kmertochar[v];

            ctx = ((ctx << 4) | uv) & E->ctx_mask;

            i = j + 1;
        }
    }
    else {
        for (i = 0; i < n - 1 && i / 2 < prefix_len;) {
            uv = cond_dist16_decode(E->ac, &E->cs0[i/2], ctx);
            u = uv >> 2;
            v = uv & 0x3;
            x->seq.s[i++] = kmertochar[u];
            x->seq.s[i++] = kmertochar[v];
            ctx = ((ctx << 4) | uv) & E->ctx_mask;
        }

        while (i < n - 1) {
            uv = cond_dist16_decode(E->ac, &E->cs, ctx);
            u = uv >> 2;
            v = uv & 0x3;
            x->seq.s[i++] = kmertochar[u];
            x->seq.s[i++] = kmertochar[v];
            ctx = ((ctx << 4) | uv) & E->ctx_mask;
        }
    }

    if (i < n && (x->qual.n == 0 || qs[i] != 0)) {
        uv = cond_dist16_decode(E->ac, &E->cs, ctx);
        v = uv & 0x3;
        x->seq.s[i] = kmertochar[v];
    }

    x->seq.s[n] = '\0';
}


static void seqenc_decode_alignment(seqenc_t* E, seq_t* x, size_t n)
{
    static size_t N = 0;
    ++N;

    while (n >= x->seq.size) fastq_expand_str(&x->seq);


    uint32_t contig_idx, spos;
    uint32_t bytes;
    uint32_t b;
    uint32_t z;

    /* decode contig index */
    bytes = 1 + dist4_decode(E->ac, &E->d_contig_idx_bytes);
    z = 0;
    for (b = 0; b < bytes; ++b) {
        z |= cond_dist256_decode(E->ac, &E->d_contig_idx, bytes - b - 1) << (8 * b);
    }
    contig_idx = z;

    assert(contig_idx < E->contig_count);

    const twobit_t* contig = E->contigs[contig_idx];
    size_t contig_len = twobit_len(contig);


    /* decode strand */
    uint8_t strand = dist2_decode(E->ac, &E->ss);


    /* decode offset */
    bytes = 1 + dist4_decode(E->ac, &E->d_contig_off_bytes);
    z = 0;
    for (b = 0; b < bytes; ++b) {
        z |= cond_dist256_decode(E->ac, &E->d_contig_off, bytes - b - 1) << (8 * b);
    }
    spos = z;


    /* decoded edit operations */
    edit_op_t op, last_op = EDIT_MATCH;

    kmer_t u, v, ctx = 0; /* nucleotide context */

    size_t i; /* position within the read */
    size_t j; /* position within the contig */

    size_t contig_pos;

    for (i = 0, j = spos; i < n;) {
        op = cond_dist4_decode(E->ac, &E->d_edit_op, last_op);

        switch (op) {
            case EDIT_MATCH:
                if (strand) u = kmer_comp1(twobit_get(contig, contig_len - j - 1));
                else        u = twobit_get(contig, j);
                x->seq.s[i] = kmertochar[u];

                if (strand) contig_pos = contig_len - j - 1;
                else        contig_pos = j;

                v = strand ? kmer_comp1(u) : u;
                E->mismatch_tally[contig_idx][4 * contig_pos + v]++;

                ctx = ((ctx << 2) | u) & E->ins_ctx_mask;
                ++i;
                ++j;
                break;

            case EDIT_MISMATCH:
                u = cond_dist4_decode(E->ac, &E->d_ins_nuc, ctx);
                x->seq.s[i] = kmertochar[u];

                if (strand) contig_pos = contig_len - j - 1;
                else        contig_pos = j;

                v = strand ? kmer_comp1(u) : u;
                E->mismatch_tally[contig_idx][4 * contig_pos + v]++;

                ctx = ((ctx << 2) | u) & E->ins_ctx_mask;
                ++i;
                ++j;
                break;

            case EDIT_Q_GAP:
                ++j;
                break;

            case EDIT_S_GAP:
                u = cond_dist4_decode(E->ac, &E->d_ins_nuc, ctx);
                x->seq.s[i] = kmertochar[u];

                ctx = ((ctx << 2) | u) & E->ins_ctx_mask;
                ++i;
                break;
        }

        last_op = ((last_op << 2) | op) & E->edit_op_mask;
    }

    x->seq.s[n] = '\0';
}


static void seqenc_decode_contigs(seqenc_t* E)
{
    seq_t* seq = fastq_alloc_seq();

    uint32_t type;
    size_t len;

    while (E->contig_count < E->expected_contigs) {
        type = dist2_decode(E->ac, &E->ms);

        /* contigs aligned to prior contigs is not currently supported. */
        assert(type == 0);

        len = E->contig_lens[E->contig_count];

        seqenc_decode_seq(E, seq, len);
        E->contigs[E->contig_count] = twobit_alloc_n(len);
        twobit_copy_n(E->contigs[E->contig_count], seq->seq.s, len);

        E->contig_count++;
    }

    fastq_free_seq(seq);
}



void seqenc_decode(seqenc_t* E, seq_t* x, size_t n)
{
    if (E->contig_count < E->expected_contigs) {
        seqenc_decode_contigs(E);
    }

    uint32_t type = dist2_decode(E->ac, &E->ms);

    if (type == 0) seqenc_decode_seq(E, x, n);
    else           seqenc_decode_alignment(E, x, n);
}



void seqenc_flush(seqenc_t* E)
{
    ac_flush_encoder(E->ac);
}


static void patch_mismatches(seqenc_t* E)
{
    size_t len;
    size_t i, j, k;

    kmer_t u;
    uint16_t max_k;

    size_t flip_cnt = 0;

    for (i = 0; i < E->contig_count; ++i) {
        len = twobit_len(E->contigs[i]);
        for (j = 0; j < len; ++j) {
            if (* (uint64_t*) (E->mismatch_tally[i] + 4 * j) > 0) {

                for (k = 0, max_k = 0; k < 4; ++k) {
                    if (E->mismatch_tally[i][4 * j + k] > 
                        E->mismatch_tally[i][4 * j + max_k])
                    {
                        max_k = k;
                    }
                }

                u = twobit_get(E->contigs[i], j);

                if (E->mismatch_tally[i][4 * j + max_k] > mismatch_patch_cutoff &&
                    E->mismatch_tally[i][4 * j + max_k] > 
                    mismatch_patch_factor * E->mismatch_tally[i][4 * j + u])
                {
                    twobit_set(E->contigs[i], j, max_k);
                    ++flip_cnt;
                }
            }
        }

        memset(E->mismatch_tally[i], 0, 4 * len * sizeof(uint16_t));
    }

    if (verbose) fprintf(stderr, "\t%zu mismatches flipped.\n", flip_cnt);
}

void seqenc_reset_decoder(seqenc_t* E)
{
    patch_mismatches(E);
    ac_reset_decoder(E->ac);
}


