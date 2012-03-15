
#include "seqenc.h"
#include "seqenc_prior.h"
#include "misc.h"
#include "ac.h"
#include "dist.h"
#include <assert.h>
#include <string.h>


/* Order of the markov chain assigning probabilities to dinucleotides. */
static const size_t k = 11;

/* Use a seperate model for the first n dinucleotides. This is primarily to
 * account for positional sequence bias that is sommon in short read sequencing.  */
static const size_t prefix_len = 4;

/* The rate at which the nucleotide markov chain is updated. */
static const size_t seq_update_rate   = 1;
static const size_t motif_update_rate = 4;

/* Initial pseudocount biasing contig motifs towards the consensus sequence.  */
static const uint16_t contig_motif_prior = 50;

enum {
    SEQENC_TYPE_SEQUENCE = 0,
    SEQENC_TYPE_ALIGNMENT
};


struct seqenc_t_
{
    /* coder */
    ac_t* ac;

    /* bitmask for two bit encoded context */
    uint32_t ctx_mask;
    uint32_t ctx0_mask;

    /* nucleotide probability given the last k nucleotides */
    cond_dist16_t cs;

    /* special case models for the first 2 * prefix_les positions. */
    cond_dist16_t cs0[5];

    /* binary distribution of unique (0) / match (1) */
    dist2_t d_type;

    /* distribution over match strand */
    dist2_t d_aln_strand;

    /* distribution over the offset into the contig */
    cond_dist256_t d_contig_off;

    /* contig motifs used to compress nucleotides in alignments */
    cond_dist4_t supercontig_motif;

    /* length of the supercontig to be decoded */
    size_t expected_supercontig_len;
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

    dist2_init(&E->d_type);
    dist2_init(&E->d_aln_strand);

    cond_dist256_init(&E->d_contig_off, 9 * 256);

    E->expected_supercontig_len = 0;
    memset(&E->supercontig_motif, 0, sizeof(cond_dist4_t));
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

    cond_dist256_free(&E->d_contig_off);
    cond_dist4_free(&E->supercontig_motif);

    free(E);
}



void seqenc_encode_twobit_seq(seqenc_t* E, const twobit_t* x)
{
    dist2_encode(E->ac, &E->d_type, SEQENC_TYPE_SEQUENCE);

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
    dist2_encode(E->ac, &E->d_type, SEQENC_TYPE_SEQUENCE);

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



void seqenc_encode_alignment(
        seqenc_t* E,
        uint32_t spos, uint8_t strand,
        const twobit_t* query)
{
    size_t qlen = twobit_len(query);
    size_t slen = E->supercontig_motif.n;

    assert(spos < slen);

    dist2_encode(E->ac, &E->d_type, SEQENC_TYPE_ALIGNMENT);
    dist2_encode(E->ac, &E->d_aln_strand, strand);
    dist_encode_uint32(E->ac, &E->d_contig_off, spos);


    kmer_t u;
    size_t i;
    if (strand) {
        for (i = 0; i < qlen; ++i) {
            u = kmer_comp1(twobit_get(query, i));
            cond_dist4_encode(E->ac, &E->supercontig_motif, slen - (spos + i) - 1, u);
        }
    }
    else {
        for (i = 0; i < qlen; ++i) {
            u = twobit_get(query, i);
            cond_dist4_encode(E->ac, &E->supercontig_motif, spos + i, u);
        }
    }
}


void seqenc_set_supercontig(seqenc_t* E, const twobit_t* supercontig)
{
    size_t len = twobit_len(supercontig);
    cond_dist4_init(&E->supercontig_motif, len);
    cond_dist4_set_update_rate(&E->supercontig_motif, motif_update_rate);
    size_t i;
    kmer_t u;
    for (i = 0; i < len; ++i) {
        u = twobit_get(supercontig, i);
        E->supercontig_motif.xss[i].xs[u].count = contig_motif_prior;
        dist4_update(&E->supercontig_motif.xss[i]);
   }
}


void seqenc_get_supercontig_consensus(seqenc_t* E, twobit_t* supercontig)
{
    size_t len = twobit_len(supercontig);
    kmer_t u;
    size_t i;
    for (i = 0; i < len; ++i) {
        u = 0;
        if (E->supercontig_motif.xss[i].xs[1].count >
            E->supercontig_motif.xss[i].xs[u].count) u = 1;
        if (E->supercontig_motif.xss[i].xs[2].count >
            E->supercontig_motif.xss[i].xs[u].count) u = 2;
        if (E->supercontig_motif.xss[i].xs[3].count >
            E->supercontig_motif.xss[i].xs[u].count) u = 3;

        twobit_set(supercontig, i, u);
    }
}


void seqenc_prepare_decoder(seqenc_t* E, uint32_t supercontig_len)
{
    E->expected_supercontig_len = supercontig_len;
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
    x->seq.n = n;
}


static void seqenc_decode_alignment(seqenc_t* E, seq_t* x, size_t qlen)
{
    while (qlen >= x->seq.size) fastq_expand_str(&x->seq);

    uint8_t  strand = dist2_decode(E->ac, &E->d_aln_strand);
    uint32_t spos   = dist_decode_uint32(E->ac, &E->d_contig_off);

    size_t slen = E->supercontig_motif.n;
    
    assert(spos < slen);

    kmer_t u;
    size_t i;
    if (strand) {
        for (i = 0; i < qlen; ++i) {
            u = cond_dist4_decode(E->ac, &E->supercontig_motif, slen - (spos + i) - 1);
            x->seq.s[i] = kmertochar[kmer_comp1(u)];
        }
    }
    else {
        for (i = 0; i < qlen; ++i) {
            u = cond_dist4_decode(E->ac, &E->supercontig_motif, spos + i);
            x->seq.s[i] = kmertochar[u];
        }
    }

    x->seq.s[qlen] = '\0';
    x->seq.n = qlen;
}


static void seqenc_decode_supercontig(seqenc_t* E)
{
    seq_t* seq = fastq_alloc_seq();

#ifdef DNDEBUG
    dist2_decode(E->ac, &E->d_type);
#else
    assert(dist2_decode(E->ac, &E->d_type) == SEQENC_TYPE_SEQUENCE);
#endif

    seqenc_decode_seq(E, seq, E->expected_supercontig_len);

    twobit_t* supercontig = twobit_alloc_n(E->expected_supercontig_len);
    twobit_copy_n(supercontig, seq->seq.s, E->expected_supercontig_len);

    seqenc_set_supercontig(E, supercontig);

    twobit_free(supercontig);
    fastq_free_seq(seq);

    E->expected_supercontig_len = 0;
}



void seqenc_decode(seqenc_t* E, seq_t* x, size_t n)
{
    if (E->expected_supercontig_len > 0) seqenc_decode_supercontig(E);

    uint32_t type = dist2_decode(E->ac, &E->d_type);

    if (type == SEQENC_TYPE_SEQUENCE) seqenc_decode_seq(E, x, n);
    else                              seqenc_decode_alignment(E, x, n);
}


size_t seqenc_finish(seqenc_t* E)
{
    return ac_finish_encoder(E->ac);
}


void seqenc_flush(seqenc_t* E)
{
    ac_flush_encoder(E->ac);
}


void seqenc_start_decoder(seqenc_t* E)
{
    ac_start_decoder(E->ac);
}


void seqenc_reset_decoder(seqenc_t* E)
{
    ac_reset_decoder(E->ac);
}


