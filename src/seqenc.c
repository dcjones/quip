
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
static const size_t seq_update_rate = 2;

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
    cond_dist4_t* contig_motifs;

    /* number of contigs alignment may originate from */
    uint32_t contig_count;

    /* how many of the leading sequences are contigs */
    uint32_t expected_contigs;

    /* lengths of the expected contigs */
    uint32_t* contig_lens;

    /* cumulative contig lengths
     * We encode a position of an alignment as an offset into a
     * super-contig formed by laying the contigs end-to-end. */
    uint32_t* cum_contig_lens;
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

    E->contig_motifs    = NULL;
    E->contig_lens      = NULL;
    E->cum_contig_lens  = NULL;
    E->contig_count     = 0;
    E->expected_contigs = 0;
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

    for (i = 0; i < E->contig_count; ++i) {
        cond_dist4_free(&E->contig_motifs[i]);
    }
    free(E->contig_motifs);
    free(E->contig_lens);
    free(E->cum_contig_lens);

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
    prefetch(x);

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
        uint32_t contig_idx, uint32_t spos, uint8_t strand,
        const twobit_t* query)
{
    cond_dist4_t* motif = &E->contig_motifs[contig_idx];
    size_t qlen = twobit_len(query);
    size_t slen = motif->n;

    dist2_encode(E->ac, &E->d_type, SEQENC_TYPE_ALIGNMENT);

    dist2_encode(E->ac, &E->d_aln_strand, strand);

    dist_encode_uint32(E->ac, &E->d_contig_off,
        E->cum_contig_lens[contig_idx] + spos);


    kmer_t u;
    size_t i;
    if (strand) {
        for (i = 0; i < qlen; ++i) {
            u = kmer_comp1(twobit_get(query, i));
            cond_dist4_encode(E->ac, motif, slen - (spos + i) - 1, u);
        }
    }
    else {
        for (i = 0; i < qlen; ++i) {
            u = twobit_get(query, i);
            cond_dist4_encode(E->ac, motif, spos + i, u);
        }
    }
}


void seqenc_set_contigs(seqenc_t* E, twobit_t** contigs, size_t n)
{
    size_t i, j;
    for (i = 0; i < E->contig_count; ++i) {
        cond_dist4_free(&E->contig_motifs[i]);
    }
    free(E->contig_motifs);
    free(E->cum_contig_lens);

    E->cum_contig_lens = malloc_or_die(n * sizeof(uint32_t));

    twobit_t*     contig;
    cond_dist4_t* motif;
    size_t len;

    E->contig_count = n;
    E->contig_motifs = malloc_or_die(n * sizeof(cond_dist4_t));

    if (n > 0) E->cum_contig_lens[0] = 0;

    for (i = 0; i < E->contig_count; ++i) {
        contig = contigs[i];
        motif  = &E->contig_motifs[i];
        len    = twobit_len(contig);

        cond_dist4_init(motif, len);
        for (j = 0; j < len; ++j) {
            motif->xss[j].xs[twobit_get(contig, j)].count =
                contig_motif_prior;
            dist4_update(&motif->xss[j]);
        }
        cond_dist4_set_update_rate(motif, 4);

        if (i + 1 < E->contig_count) {
            E->cum_contig_lens[i + 1] = E->cum_contig_lens[i] + len;
        }
    }
}


void seqenc_get_contig_consensus(seqenc_t* E, twobit_t** contigs)
{
    size_t i, j;
    twobit_t*     contig;
    cond_dist4_t* motif;
    size_t len;
    kmer_t u;

    for (i = 0; i < E->contig_count; ++i) {
        contig = contigs[i];
        motif  = &E->contig_motifs[i];
        len    = twobit_len(contig);;

        for (j = 0; j < len; ++j) {
            u = 0;
            if (motif->xss[j].xs[1].count >
                motif->xss[j].xs[u].count) u = 1;
            if (motif->xss[j].xs[2].count >
                motif->xss[j].xs[u].count) u = 2;
            if (motif->xss[j].xs[3].count >
                motif->xss[j].xs[u].count) u = 3;

            twobit_set(contig, j, u);
        }
    }
}


void seqenc_prepare_decoder(seqenc_t* E, uint32_t n, const uint32_t* lens)
{
    size_t i;

    free(E->contig_lens);
    E->contig_lens = NULL;

    E->contig_count = 0;
    E->expected_contigs = n;

    E->contig_lens = malloc_or_die(n * sizeof(uint32_t));

    for (i = 0; i < n; ++i) {
        E->contig_lens[i] = lens[i];
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
    x->seq.n = n;
}


static uint32_t decode_contig_idx(seqenc_t* E, uint32_t supercontig_off)
{
    /* find the contig index by binary search */
    int i, j, k;
    i = 0;
    j = E->contig_count;

    do {
        k = (i + j) / 2;
        if (supercontig_off < E->cum_contig_lens[k]) j = k;
        else i = k;
    } while (i + 1 < j);

    assert(i == (int) E->contig_count - 1 || supercontig_off < E->cum_contig_lens[i + 1]);
    assert(supercontig_off >= E->cum_contig_lens[i]);

    return i;
}


static void seqenc_decode_alignment(seqenc_t* E, seq_t* x, size_t qlen)
{
    while (qlen >= x->seq.size) fastq_expand_str(&x->seq);

    uint32_t contig_idx, spos;
    uint32_t z;

    uint8_t strand = dist2_decode(E->ac, &E->d_aln_strand);

    /* contig index and offset */
    z = dist_decode_uint32(E->ac, &E->d_contig_off);
    contig_idx = decode_contig_idx(E, z);
    spos = z - E->cum_contig_lens[contig_idx];

    /* contig motif */
    if (contig_idx >= E->contig_count) {
        fprintf(stderr, "Error: contig index out of bounds.\n");
        exit(EXIT_FAILURE);
    }
    cond_dist4_t* motif = &E->contig_motifs[contig_idx];
    size_t slen = motif->n;

    kmer_t u;
    size_t i;
    if (strand) {
        for (i = 0; i < qlen; ++i) {
            u = cond_dist4_decode(E->ac, motif, slen - (spos + i) - 1);
            x->seq.s[i] = kmertochar[kmer_comp1(u)];
        }
    }
    else {
        for (i = 0; i < qlen; ++i) {
            u = cond_dist4_decode(E->ac, motif, spos + i);
            x->seq.s[i] = kmertochar[u];
        }
    }

    x->seq.s[qlen] = '\0';
    x->seq.n = qlen;
}


static void seqenc_decode_contigs(seqenc_t* E)
{
    seq_t* seq = fastq_alloc_seq();

    uint32_t type;
    size_t len;

    twobit_t** contigs = NULL;
    contigs = malloc_or_die(E->expected_contigs * sizeof(twobit_t*));
    memset(contigs, 0, E->expected_contigs * sizeof(twobit_t*));

    uint32_t cnt = 0;

    while (cnt < E->expected_contigs) {
        type = dist2_decode(E->ac, &E->d_type);

        /* contigs aligned to prior contigs is not currently supported. */
        assert(type == SEQENC_TYPE_SEQUENCE);

        len = E->contig_lens[cnt];

        seqenc_decode_seq(E, seq, len);
        contigs[cnt] = twobit_alloc_n(len);
        twobit_copy_n(contigs[cnt], seq->seq.s, len);

        cnt++;
    }

    fastq_free_seq(seq);

    /* initialize contig motifs */
    seqenc_set_contigs(E, contigs, cnt);

    uint32_t i;
    for (i = 0; i < cnt; ++i) {
        twobit_free(contigs[i]);
    }
    free(contigs);

    if (quip_verbose) fprintf(stderr, "read %u contigs\n", (unsigned int) cnt);
}



void seqenc_decode(seqenc_t* E, seq_t* x, size_t n)
{
    if (E->contig_count < E->expected_contigs) {
        seqenc_decode_contigs(E);
    }

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


