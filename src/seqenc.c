
#include "seqenc.h"
#include "seqenc_prior.h"
#include "misc.h"
#include "ac.h"
#include "dist.h"
#include "seqmap.h"
#include "sam/bam.h"
#include <assert.h>
#include <string.h>


/* Order of the markov chain assigning probabilities to dinucleotides. */
static const size_t k = 11;

/* Order of the markov chain assigning probabilities to the N mask. */
static const size_t k_nmask = 6;

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

enum {
    SEQENC_REF_MATCH = 0,
    SEQENC_REF_MISMATCH
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

    /* N probability, given the last k_nmask bits. */
    cond_dist2_t d_nmask;
    uint32_t nmask_ctx_mask;

    /* binary distribution of unique (0) / match (1) */
    dist2_t d_type;

    /* distribution over match strand */
    dist2_t d_aln_strand;

    /* distribution over the offset into the contig */
    cond_dist256_t d_contig_off;

    /* distribution over match/mismatces in reference alignments */
    dist2_t d_ref_match;

    /* distribution over inserted nucleotides in reference alignment */
    dist4_t d_ref_ins_nuc;

    /* contig motifs used to compress nucleotides in alignments */
    cond_dist4_t supercontig_motif;

    /* reference sequence set, or NULL if none */
    const seqmap_t* ref;

    /* distribution used to compress the "extra" fields in
     * the short_read_t structure. */

    /* flags */
    cond_dist256_t  d_ext_flags;

    /* sequence index */
    cond_dist256_t  d_ext_seqidx;

    /* genomic position, conditioned on sequence index */
    cond_dist256_t* d_ext_pos;

    /* map quality */
    dist256_t       d_ext_map_qual;

    /* cigar operation, conditioned on the previous cigar operation */
    cond_dist16_t   d_ext_cigar_op;

    /* cigar length, conditioned on the operation */
    cond_dist256_t  d_ext_cigar_len[9];

    /* whether the mate is aligned to the same sequence */
    dist2_t         d_ext_mate_sameseq;

    /* template length */
    cond_dist256_t  d_ext_tlen;

    /* first-order markov chain for aux data */
    cond_dist256_t  d_ext_aux_len;
    cond_dist256_t  d_ext_aux;
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


static void seqenc_init(seqenc_t* E, const seqmap_t* ref)
{
    E->ref = ref;

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

    N = 1 << k_nmask;
    E->nmask_ctx_mask = N -1 ;
    cond_dist2_init(&E->d_nmask, N);

    dist2_init(&E->d_type);
    dist2_init(&E->d_aln_strand);

    cond_dist256_init(&E->d_contig_off, 9 * 256);

    dist2_init(&E->d_ref_match);
    dist4_init(&E->d_ref_ins_nuc);

    memset(&E->supercontig_motif, 0, sizeof(cond_dist4_t));

    size_t refsize = ref == NULL ? 0 : seqmap_size(ref);

    cond_dist256_init(&E->d_ext_flags,  9 * 256);
    cond_dist256_init(&E->d_ext_seqidx, 9 * 256);

    E->d_ext_pos = malloc_or_die(refsize * sizeof(cond_dist256_t));
    for (i = 0; i < refsize; ++i) {
        cond_dist256_init(&E->d_ext_pos[i], 9 * 256);
    }

    dist256_init(&E->d_ext_map_qual);
    cond_dist16_init(&E->d_ext_cigar_op, 10);

    for (i = 0; i < 9; ++i) {
        cond_dist256_init(&E->d_ext_cigar_len[i], 9 * 256);
    }

    dist2_init(&E->d_ext_mate_sameseq);
    cond_dist256_init(&E->d_ext_tlen, 9 * 256);
    cond_dist256_init(&E->d_ext_aux_len, 9 * 256);
    cond_dist256_init(&E->d_ext_aux, 256);
}


seqenc_t* seqenc_alloc_encoder(quip_writer_t writer, void* writer_data, const seqmap_t* ref)
{
    seqenc_t* E = malloc_or_die(sizeof(seqenc_t));

    E->ac = ac_alloc_encoder(writer, writer_data);

    seqenc_init(E, ref);

    return E;
}


seqenc_t* seqenc_alloc_decoder(quip_reader_t reader, void* reader_data, const seqmap_t* ref)
{
    seqenc_t* E = malloc_or_die(sizeof(seqenc_t));

    E->ac = ac_alloc_decoder(reader, reader_data);

    seqenc_init(E, ref);

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

    cond_dist2_free(&E->d_nmask);

    cond_dist256_free(&E->d_contig_off);
    cond_dist4_free(&E->supercontig_motif);


    cond_dist256_free(&E->d_ext_flags);
    cond_dist256_free(&E->d_ext_seqidx);

    size_t refsize = E->ref == NULL ? 0 : seqmap_size(E->ref);
    for (i = 0; i < refsize; ++i) {
        cond_dist256_free(&E->d_ext_pos[i]);
    }
    free(E->d_ext_pos);

    cond_dist16_free(&E->d_ext_cigar_op);

    for (i = 0; i < 9; ++i) {
        cond_dist256_free(&E->d_ext_cigar_len[i]);
    }

    cond_dist256_free(&E->d_ext_tlen);
    cond_dist256_free(&E->d_ext_aux_len);
    cond_dist256_free(&E->d_ext_aux);

    free(E);
}


void seqenc_encode_extras(seqenc_t* E, const short_read_t* x)
{
    size_t i;

    dist_encode_uint32(E->ac, &E->d_ext_flags, x->flags);
    dist256_encode(E->ac, &E->d_ext_map_qual, x->map_qual);
    dist_encode_uint32(E->ac, &E->d_ext_tlen, x->tlen);

    /* aux data */
    dist_encode_uint32(E->ac, &E->d_ext_aux_len, x->aux.n);
    unsigned char last = '\0';
    for (i = 0; i < x->aux.n; ++i) {
        cond_dist256_encode(E->ac, &E->d_ext_aux, x->aux.s[i], last);
        last = x->aux.s[i];
    }

    if ((x->flags & BAM_FUNMAP) == 0) {
        if (E->ref) {
            /* TODO: lookup sequence index in reference */
        }
        else {
            /* TODO: lookup sequence index in hash table */
        }
        /* TODO: encode sequence index */

        /* TODO: encode position */

        uint8_t last_op = 9;
        for (i = 0; i < x->cigar.n; ++i) {
            cond_dist16_encode(E->ac, &E->d_ext_cigar_op, x->cigar.ops[i], last_op);
            dist_encode_uint32(E->ac, &E->d_ext_cigar_len[x->cigar.ops[i]], x->cigar.lens[i]);
        }
    }

    if ((x->flags & BAM_FMUNMAP) == 0) {
        if (E->ref) {
            /* TODO: lookup sequence index in reference */
        }
        else {
            /* TODO: lookup sequence index in hash table */
        }
        /* TODO: encode whether it matches the main sequence */
        /*       if not, encode the sequence index */

        /* TODO: encode position */
    }
}


void seqenc_encode_twobit_seq(seqenc_t* E, const twobit_t* x)
{
    dist2_encode(E->ac, &E->d_type, SEQENC_TYPE_SEQUENCE);

    size_t n = twobit_len(x);
    if (n == 0) return;

   
    kmer_t uv;
    uint32_t ctx = 0;
    size_t i;

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

    /* encode N mask */
    for (i = 0; i < n; ++i) {
        cond_dist2_encode(E->ac, &E->d_nmask, 0, 0);
    }
}



void seqenc_encode_char_seq(seqenc_t* E, const uint8_t* x, size_t len)
{
    dist2_encode(E->ac, &E->d_type, SEQENC_TYPE_SEQUENCE);

    kmer_t uv;
    uint32_t ctx = 0;
    size_t i;

    /* encode leading positions. */
    for (i = 0; i < len - 1 && i / 2 < prefix_len; i += 2) {
        uv = (chartokmer[x[i]] << 2) | chartokmer[x[i + 1]];
        cond_dist16_encode(E->ac, &E->cs0[i/2], ctx, uv);
        ctx = ((ctx << 4) | uv) & E->ctx_mask;
    }

    /* encode trailing positions. */
    for (; i < len - 1; i += 2) {
        uv = (chartokmer[x[i]] << 2) | chartokmer[x[i + 1]];
        cond_dist16_encode(E->ac, &E->cs, ctx, uv);
        ctx = ((ctx << 4) | uv) & E->ctx_mask;
    }

    /* handle odd read lengths */
    if (i == len - 1) {
        uv = chartokmer[x[i]];
        cond_dist16_encode(E->ac, &E->cs, ctx, uv);
    }

    /* encode N mask */
    uint32_t nmask_ctx = 0;
    for (i = 0; i < len; ++i) {
        if (x[i] == 'N') {
            cond_dist2_encode(E->ac, &E->d_nmask, nmask_ctx, 1);
            nmask_ctx = ((nmask_ctx << 1) | 1) & E->nmask_ctx_mask;
        }
        else {
            cond_dist2_encode(E->ac, &E->d_nmask, nmask_ctx, 0);
            nmask_ctx = (nmask_ctx << 1) & E->nmask_ctx_mask;
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


void seqenc_encode_reference_alignment(
        seqenc_t* E, const short_read_t* r)
{
    size_t refseq_idx;
    const twobit_t* refseq = seqmap_get(E->ref, (const char*) r->seqname.s, &refseq_idx);

    if (refseq == NULL) {
        fprintf(stderr,
            "A read was aligned to sequence %s, which was not found in the reference.\n",
            r->seqname.s);
        exit(EXIT_FAILURE);
    }

    uint32_t ref_pos  = r->pos;
    uint32_t read_pos = 0;

    size_t i = 0; /* cigar operation */
    size_t j; /* position within the cigar op */

    kmer_t x; /* read nucleotide */
    kmer_t y; /* reference nucleotide */

    for (i = 0; i < r->cigar.n; ++i) {
        switch (r->cigar.ops[i]) {
            case BAM_CEQUAL:
            case BAM_CDIFF:
            case BAM_CMATCH:
                for (j = 0; j < r->cigar.lens[i]; ++j) {
                    if (r->strand) {
                        x = kmer_comp1(chartokmer[r->seq.s[r->seq.n - 1 - read_pos]]);
                    }
                    else x = chartokmer[r->seq.s[read_pos]];
                    
                    // x = chartokmer[r->seq.s[read_pos]];
                    y = twobit_get(refseq, ref_pos);

                    if (x == y) {
                        dist2_encode(E->ac, &E->d_ref_match, SEQENC_REF_MATCH);
                    }
                    else {
                        dist2_encode(E->ac, &E->d_ref_match, SEQENC_REF_MISMATCH);
                        dist4_encode(E->ac, &E->d_ref_ins_nuc, chartokmer[r->seq.s[read_pos]]);
                    }

                    read_pos++;
                    ref_pos++;
                }
                break;

            case BAM_CINS:
                for (j = 0; j < r->cigar.lens[j]; ++j) {
                    dist4_encode(E->ac, &E->d_ref_ins_nuc, chartokmer[r->seq.s[read_pos]]);
                    read_pos++;
                }
                break;

            case BAM_CDEL:
                ref_pos += r->cigar.lens[j];
                break;

            case BAM_CREF_SKIP:
                ref_pos += r->cigar.lens[j];
                break;

            case BAM_CSOFT_CLIP:
                for (j = 0; j < r->cigar.lens[j]; ++j) {
                    dist4_encode(E->ac, &E->d_ref_ins_nuc, chartokmer[r->seq.s[read_pos]]);
                    read_pos++;
                    ref_pos++;
                }
                break;

            case BAM_CHARD_CLIP:
                read_pos += r->cigar.lens[j];
                ref_pos  += r->cigar.lens[j];
                break;

            case BAM_CPAD:
                /* TODO: I have no fucking clue what this is supposed to mean. */
                fprintf(stderr, "Unsupported cigar operation.\n");
                exit(EXIT_FAILURE);
                break;
        }
    }

    if (read_pos != r->seq.n) {
        fprintf(stderr, "Cigar operations do not account for full read length.\n");
        exit(EXIT_FAILURE);
    }

    /* encode N mask */
    uint32_t nmask_ctx = 0;
    for (i = 0; i < r->seq.n; ++i) {
        if (r->seq.s[i] == 'N') {
            cond_dist2_encode(E->ac, &E->d_nmask, nmask_ctx, 1);
            nmask_ctx = ((nmask_ctx << 1) | 1) & E->nmask_ctx_mask;
        }
        else {
            cond_dist2_encode(E->ac, &E->d_nmask, nmask_ctx, 0);
            nmask_ctx = (nmask_ctx << 1) & E->nmask_ctx_mask;
        }
    }
}


void seqenc_set_supercontig(seqenc_t* E, const twobit_t* supercontig)
{
    size_t len = twobit_len(supercontig);
    if (len == 0) return;

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


static void seqenc_decode_seq(seqenc_t* E, short_read_t* x, size_t n)
{
    if (n == 0) return;
    str_reserve(&x->seq, n + 1);

    kmer_t uv, u, v;
    uint32_t ctx = 0;
    size_t i;

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

    if (i == n - 1) {
        uv = cond_dist16_decode(E->ac, &E->cs, ctx);
        v = uv & 0x3;
        x->seq.s[i] = kmertochar[v];
    }

    uint32_t nmask_ctx = 0;
    for (i = 0; i < n; ++i) {
        nmask_ctx = (nmask_ctx << 1) | cond_dist2_decode(E->ac, &E->d_nmask, nmask_ctx);
        nmask_ctx &= E->nmask_ctx_mask;

        if (nmask_ctx & 0x1) x->seq.s[i] = 'N';
    }


    x->seq.s[n] = '\0';
    x->seq.n = n;
}


static void seqenc_decode_alignment(seqenc_t* E, short_read_t* x, size_t qlen)
{
    str_reserve(&x->seq, qlen + 1);

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



void seqenc_decode(seqenc_t* E, short_read_t* x, size_t n)
{
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


