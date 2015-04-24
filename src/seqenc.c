
#include "seqenc.h"
#include "misc.h"
#include "ac.h"
#include "dist.h"
#include "seqmap.h"
#include "strmap.h"
#include "sam/bam.h"
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

    dist2_t* d_nmask;
    size_t nmask_n; /* maximum read length supported by d_nmask */

    /* binary distribution of unique (0) / match (1) */
    dist2_t d_type;

    /* distribution over match strand */
    dist2_t d_aln_strand;

    /* distribution over the offset into the contig */
    uint32_enc_t d_contig_off;

    /* distribution over match/mismatces in reference alignments */
    dist2_t d_ref_match;

    /* distribution over inserted nucleotides in reference alignment */
    dist4_t d_ref_ins_nuc;

    /* contig motifs used to compress nucleotides in alignments */
    cond_dist4_t supercontig_motif;

    /* reference sequence set, or NULL if none */
    const seqmap_t* ref;

    /* temporary space to compute reverse complements */
    str_t tmpseq;

    /* assign sequential indexes to sequences */
    strmap_t* seq_index;

    /* distribution used to compress the "extra" fields in
     * the short_read_t structure. */

    /* flags */
    uint32_enc_t d_ext_flags;

    /* sequence name */
    cond_dist128_t  d_ext_seqname;

    /* genomic position, conditioned on sequence index */
    uint32_enc_t*   d_ext_pos;

    /* genomic position as offset from the previous */
    uint32_t        last_ref_pos;
    dist2_t         d_ext_pos_off_flag;
    dist256_t       d_ext_pos_off;

    /* map quality */
    dist256_t       d_ext_map_qual;

    /* number of cigar operations */
    uint32_enc_t    d_ext_cigar_n;

    /* cigar operation, conditioned on the previous cigar operation */
    cond_dist16_t   d_ext_cigar_op;

    /* cigar length, conditioned on the operation */
    uint32_enc_t    d_ext_cigar_len[9];

    /* whether the mate is aligned to the same sequence */
    dist2_t         d_ext_mate_sameseq;

    /* template length */
    uint32_enc_t    d_ext_tlen;
};


static void seqenc_init(seqenc_t* E, const seqmap_t* ref)
{
    E->ref = ref;
    str_init(&E->tmpseq);

    size_t N = 1 << (2 * k);
    E->ctx_mask = N - 1;

    cond_dist16_init(&E->cs, N);
    cond_dist16_set_update_rate(&E->cs, seq_update_rate);

    size_t i;
    for (i = 0; i < prefix_len; ++i) {
        cond_dist16_init(&E->cs0[i], 1 << (4 * i));
        cond_dist16_set_update_rate(&E->cs0[i], seq_update_rate);
    }

    E->d_nmask = NULL;
    E->nmask_n = 0;

    dist2_init(&E->d_type);
    dist2_init(&E->d_aln_strand);

    uint32_enc_init(&E->d_contig_off);

    dist2_init(&E->d_ref_match);
    dist4_init(&E->d_ref_ins_nuc);

    memset(&E->supercontig_motif, 0, sizeof(cond_dist4_t));

    uint32_enc_init(&E->d_ext_flags);
    cond_dist128_init(&E->d_ext_seqname, 128);

    E->d_ext_pos = NULL;

    E->last_ref_pos = 0;
    dist2_init(&E->d_ext_pos_off_flag);
    dist256_init(&E->d_ext_pos_off);

    dist256_init(&E->d_ext_map_qual);
    uint32_enc_init(&E->d_ext_cigar_n);
    cond_dist16_init(&E->d_ext_cigar_op, 10);

    for (i = 0; i < 9; ++i) {
        uint32_enc_init(&E->d_ext_cigar_len[i]);
    }

    dist2_init(&E->d_ext_mate_sameseq);
    uint32_enc_init(&E->d_ext_tlen);

    E->seq_index = strmap_alloc();
}


static void reserve_nmask(seqenc_t* E, size_t readlen)
{
    if (readlen > E->nmask_n) {
        E->d_nmask = realloc_or_die(E->d_nmask, readlen * sizeof(dist2_t));

        size_t i;
        for (i = E->nmask_n; i < readlen; ++i) {
            dist2_init(&E->d_nmask[i]);
        }

        E->nmask_n = readlen;
    }
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

    str_free(&E->tmpseq);

    ac_free(E->ac);
    cond_dist16_free(&E->cs);

    size_t i;
    for (i = 0; i < prefix_len; ++i) {
        cond_dist16_free(&E->cs0[i]);
    }

    free(E->d_nmask);

    uint32_enc_free(&E->d_contig_off);
    cond_dist4_free(&E->supercontig_motif);

    uint32_enc_free(&E->d_ext_flags);
    cond_dist128_free(&E->d_ext_seqname);

    size_t refsize = strmap_size(E->seq_index);
    for (i = 0; i < refsize; ++i) {
        cond_dist256_free(&E->d_ext_pos[i]);
    }
    free(E->d_ext_pos);

    cond_dist16_free(&E->d_ext_cigar_op);

    uint32_enc_free(&E->d_ext_cigar_n);
    for (i = 0; i < 9; ++i) {
        uint32_enc_free(&E->d_ext_cigar_len[i]);
    }

    uint32_enc_free(&E->d_ext_tlen);

    strmap_free(E->seq_index);

    free(E);
}


static uint32_t get_seq_idx(seqenc_t* E, const str_t* seqname)
{
    uint32_t n   = strmap_size(E->seq_index);
    uint32_t idx = strmap_get(E->seq_index, seqname);

    if (idx >= n) {
        E->d_ext_pos = realloc_or_die(E->d_ext_pos, (n + 1) * sizeof(cond_dist256_t));
        cond_dist256_init(&E->d_ext_pos[n], 9 * 256);
    }

    return idx;
}


static void encode_seqname(seqenc_t* E, const str_t* seqname)
{
    unsigned char last = '\0';
    size_t i;
    for (i = 0; i < seqname->n; ++i) {
        cond_dist128_encode(E->ac, &E->d_ext_seqname, last, seqname->s[i]);
        last = seqname->s[i];
    }
    cond_dist128_encode(E->ac, &E->d_ext_seqname, last, '\0');
}


static void decode_seqname(seqenc_t* E, str_t* seqname)
{
    seqname->n = 0;
    unsigned char last = '\0';
    size_t i = 0;
    do {
        str_reserve_extra(seqname, 2);
        seqname->s[i] = \
            cond_dist128_decode(E->ac, &E->d_ext_seqname, last);
        last = seqname->s[i];
        ++seqname->n;
        ++i;

    } while (last != '\0');
    seqname->n = i - 1;
}


void seqenc_encode_extras(seqenc_t* E, const short_read_t* x, uint8_t quip_version)
{
    uint32_enc_encode(E->ac, &E->d_ext_flags, x->flags);
    dist256_encode(E->ac, &E->d_ext_map_qual, x->map_qual);
    uint32_enc_encode(E->ac, &E->d_ext_tlen, x->tlen);

    uint32_t seqidx = 0;
    if ((x->flags & BAM_FUNMAP) == 0) {
        encode_seqname(E, &x->seqname);
        seqidx = get_seq_idx(E, &x->seqname);

        if (x->pos < E->last_ref_pos || x->pos - E->last_ref_pos >= 256) {
            dist2_encode(E->ac, &E->d_ext_pos_off_flag, 0);
            uint32_enc_encode(E->ac, &E->d_ext_pos[seqidx], x->pos);
        }
        else {
            dist2_encode(E->ac, &E->d_ext_pos_off_flag, 1);
            dist256_encode(E->ac, &E->d_ext_pos_off, x->pos - E->last_ref_pos);
        }

        E->last_ref_pos = x->pos;

        uint32_t cigarlen = 0;
        uint8_t last_op = 9;
        size_t i;

        uint32_enc_encode(E->ac, &E->d_ext_cigar_n, x->cigar.n);
        for (i = 0; i < x->cigar.n; ++i) {
            cond_dist16_encode(E->ac, &E->d_ext_cigar_op, last_op, x->cigar.ops[i]);
            uint32_enc_encode(E->ac, &E->d_ext_cigar_len[x->cigar.ops[i]], x->cigar.lens[i]);
            last_op = x->cigar.ops[i];

            if (x->cigar.ops[i] != BAM_CDEL &&
                x->cigar.ops[i] != BAM_CREF_SKIP &&
                x->cigar.ops[i] != BAM_CHARD_CLIP)
            {
                cigarlen += x->cigar.lens[i];
            }
        }

        if (cigarlen != x->seq.n) {
            quip_error("Cigar operations do not account for full read length.");
        }
    }

    if ((x->flags & BAM_FMUNMAP) == 0) {
        if (quip_version >= 4) {
            if ((x->flags & BAM_FUNMAP) == 0) {
                if (strcmp((char*) x->seqname.s, (char*) x->mate_seqname.s) == 0) {
                    dist2_encode(E->ac, &E->d_ext_mate_sameseq, 1);
                }
                else {
                    dist2_encode(E->ac, &E->d_ext_mate_sameseq, 0);
                    encode_seqname(E, &x->mate_seqname);
                }
            }
            else {
                encode_seqname(E, &x->mate_seqname);
            }
            seqidx = get_seq_idx(E, &x->mate_seqname);
        }
        else {
            if ((x->mate_seqname.n == 1 && x->mate_seqname.s[0] == '=') ||
                ((x->flags & BAM_FUNMAP) == 0 && strcmp((char*) x->seqname.s, (char*) x->mate_seqname.s) == 0))
            {
                dist2_encode(E->ac, &E->d_ext_mate_sameseq, 1);
            }
            else {
                dist2_encode(E->ac, &E->d_ext_mate_sameseq, 0);
                encode_seqname(E, &x->mate_seqname);
                seqidx = get_seq_idx(E, &x->mate_seqname);
            }
        }

        uint32_enc_encode(E->ac, &E->d_ext_pos[seqidx], x->mate_pos);
    }
}


void seqenc_decode_extras(seqenc_t* E, short_read_t* x, size_t seqlen,
                          uint8_t quip_version)
{
    x->flags    = uint32_enc_decode(E->ac, &E->d_ext_flags);
    x->strand   = (x->flags & BAM_FREVERSE) ? 1 : 0;
    x->map_qual = dist256_decode(E->ac, &E->d_ext_map_qual);
    x->tlen     = uint32_enc_decode(E->ac, &E->d_ext_tlen);

    x->cigar.n = 0;
    uint32_t seqidx = 0;
    if ((x->flags & BAM_FUNMAP) == 0) {
        decode_seqname(E, &x->seqname);
        seqidx = get_seq_idx(E, &x->seqname);

        if (dist2_decode(E->ac, &E->d_ext_pos_off_flag)) {
            x->pos =
                E->last_ref_pos + dist256_decode(E->ac, &E->d_ext_pos_off);
        }
        else {
            x->pos = uint32_enc_decode(E->ac, &E->d_ext_pos[seqidx]);
        }

        E->last_ref_pos = x->pos;

        uint8_t last_op = 9;
        size_t i = 0;
        uint32_t cigarlen = 0;
        x->cigar.n = uint32_enc_decode(E->ac, &E->d_ext_cigar_n);
        cigar_reserve(&x->cigar, x->cigar.n);

        for (i = 0; i < x->cigar.n; ++i) {
            x->cigar.ops[i] = cond_dist16_decode(E->ac, &E->d_ext_cigar_op, last_op);
            x->cigar.lens[i] = uint32_enc_decode(E->ac, &E->d_ext_cigar_len[x->cigar.ops[i]]);

            if (x->cigar.ops[i] != BAM_CDEL &&
                x->cigar.ops[i] != BAM_CREF_SKIP &&
                x->cigar.ops[i] != BAM_CHARD_CLIP)
            {
                cigarlen += x->cigar.lens[i];
            }

            last_op = x->cigar.ops[i];
        }

        if (cigarlen != seqlen) {
            quip_error("Cigar operations do not account for full read length.");
        }
    }

    if ((x->flags & BAM_FMUNMAP) == 0) {
        if (quip_version >= 4) {
            if ((x->flags & BAM_FUNMAP) == 0) {
                if (dist2_decode(E->ac, &E->d_ext_mate_sameseq)) {
                    str_copy(&x->mate_seqname, &x->seqname);
                }
                else {
                    decode_seqname(E, &x->mate_seqname);
                }
            }
            else {
                decode_seqname(E, &x->mate_seqname);
            }

            seqidx = get_seq_idx(E, &x->mate_seqname);
        }
        else {
            if (dist2_decode(E->ac, &E->d_ext_mate_sameseq)) {
                str_copy(&x->mate_seqname, &x->seqname);
            }
            else {
                decode_seqname(E, &x->mate_seqname);
                seqidx = get_seq_idx(E, &x->mate_seqname);
            }
        }
        x->mate_pos = uint32_enc_decode(E->ac, &E->d_ext_pos[seqidx]);
    }
}


void seqenc_encode_twobit_seq(seqenc_t* E, const unsigned char* x_str, const twobit_t* x)
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
    reserve_nmask(E, n);
    for (i = 0; i < n; ++i) {
        dist2_encode(E->ac, &E->d_nmask[i],
            x_str[i] == 'N' ? 1 : 0);
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
    reserve_nmask(E, len);
    for (i = 0; i < len; ++i) {
        dist2_encode(E->ac, &E->d_nmask[i],
            x[i] == 'N' ? 1 : 0);
    }

}



void seqenc_encode_alignment(
        seqenc_t* E,
        uint32_t spos, uint8_t strand,
        const unsigned char* query_str,
        const twobit_t* query)
{
    size_t qlen = twobit_len(query);
    size_t slen = E->supercontig_motif.n;

    assert(spos < slen);

    dist2_encode(E->ac, &E->d_type, SEQENC_TYPE_ALIGNMENT);

    /* encode N mask */
    size_t i;
    reserve_nmask(E, qlen);
    for (i = 0; i < qlen; ++i) {
        dist2_encode(E->ac, &E->d_nmask[i],
            query_str[i] == 'N' ? 1 : 0);
    }

    dist2_encode(E->ac, &E->d_aln_strand, strand);
    uint32_enc_encode(E->ac, &E->d_contig_off, spos);

    kmer_t u;
    if (strand) {
        for (i = 0; i < qlen; ++i) {
            if (query_str[i] == 'N') continue;
            u = kmer_comp1(twobit_get(query, i));
            cond_dist4_encode(E->ac, &E->supercontig_motif, slen - (spos + i) - 1, u);
        }
    }
    else {
        for (i = 0; i < qlen; ++i) {
            if (query_str[i] == 'N') continue;
            u = twobit_get(query, i);
            cond_dist4_encode(E->ac, &E->supercontig_motif, spos + i, u);
        }
    }
}


void seqenc_encode_reference_alignment(
        seqenc_t* E, const short_read_t* r)
{
    const twobit_t* refseq = seqmap_get(E->ref, (const char*) r->seqname.s);

    if (refseq == NULL) {
        quip_error(
            "A read was aligned to sequence %s, which was not found in the reference.\n",
            r->seqname.s);
    }

    str_copy(&E->tmpseq, &r->seq);
    if (r->strand) {
        str_revcomp(E->tmpseq.s, E->tmpseq.n);
    }

    /* encode N mask */
    size_t i;
    reserve_nmask(E, r->seq.n);
    for (i = 0; i < r->seq.n; ++i) {
        dist2_encode(E->ac, &E->d_nmask[i],
            E->tmpseq.s[i] == 'N' ? 1 : 0);
    }

    uint32_t ref_pos   = r->pos;
    uint32_t read_pos  = 0;

    i = 0;    /* cigar operation */
    size_t j; /* position within the cigar op */

    kmer_t x; /* read nucleotide */
    kmer_t y; /* reference nucleotide */

    for (i = 0; i < r->cigar.n; ++i) {
        switch (r->cigar.ops[i]) {
            case BAM_CEQUAL:
            case BAM_CDIFF:
            case BAM_CMATCH:
                for (j = 0; j < r->cigar.lens[i]; ++j, ++read_pos, ++ref_pos) {
                    if (E->tmpseq.s[read_pos] == 'N') continue;
                       
                    x = chartokmer[E->tmpseq.s[read_pos]];
                    y = twobit_get(refseq, ref_pos);

                    if (x == y) {
                        dist2_encode(E->ac, &E->d_ref_match, SEQENC_REF_MATCH);
                    }
                    else {
                        dist2_encode(E->ac, &E->d_ref_match, SEQENC_REF_MISMATCH);
                        dist4_encode(E->ac, &E->d_ref_ins_nuc, chartokmer[E->tmpseq.s[read_pos]]);
                    }
                }
                break;

            case BAM_CINS:
                for (j = 0; j < r->cigar.lens[i]; ++j, ++read_pos) {
                    if (E->tmpseq.s[read_pos] == 'N') continue;
                    dist4_encode(E->ac, &E->d_ref_ins_nuc, chartokmer[E->tmpseq.s[read_pos]]);
                }
                break;

            case BAM_CDEL:
                ref_pos += r->cigar.lens[i];
                break;

            case BAM_CREF_SKIP:
                ref_pos += r->cigar.lens[i];
                break;

            case BAM_CSOFT_CLIP:
                for (j = 0; j < r->cigar.lens[i]; ++j, ++read_pos) {
                    if (E->tmpseq.s[read_pos] == 'N') continue;
                    dist4_encode(E->ac, &E->d_ref_ins_nuc, chartokmer[E->tmpseq.s[read_pos]]);
                }
                break;

            case BAM_CHARD_CLIP:
                ref_pos += r->cigar.lens[i];
                break;

            case BAM_CPAD:
                quip_error("Cigar PAD operation is unsupported.");
                break;
        }
    }

    if (read_pos != r->seq.n) {
        quip_error("Cigar operations do not account for full read length.");
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

    /* decode N mask */
    reserve_nmask(E, n);
    for (i = 0; i < n; ++i) {
        if (dist2_decode(E->ac, &E->d_nmask[i])) x->seq.s[i] = 'N';
    }

    x->seq.s[n] = '\0';
    x->seq.n = n;
}



static void seqenc_decode_alignment(seqenc_t* E, short_read_t* x, size_t qlen)
{
    str_reserve(&x->seq, qlen + 1);
    memset(x->seq.s, '\0', qlen + 1);

    /* decode N mask */
    size_t i;
    reserve_nmask(E, qlen);
    for (i = 0; i < qlen; ++i) {
        if (dist2_decode(E->ac, &E->d_nmask[i])) x->seq.s[i] = 'N';
    }

    uint8_t  strand = dist2_decode(E->ac, &E->d_aln_strand);
    uint32_t spos   = uint32_enc_decode(E->ac, &E->d_contig_off);
    size_t slen = E->supercontig_motif.n;
   
    assert(spos < slen);

    kmer_t u;
    if (strand) {
        for (i = 0; i < qlen; ++i) {
            if (x->seq.s[i] == 'N') continue;
            u = cond_dist4_decode(E->ac, &E->supercontig_motif, slen - (spos + i) - 1);
            x->seq.s[i] = kmertochar[kmer_comp1(u)];
        }
    }
    else {
        for (i = 0; i < qlen; ++i) {
            if (x->seq.s[i] == 'N') continue;
            u = cond_dist4_decode(E->ac, &E->supercontig_motif, spos + i);
            x->seq.s[i] = kmertochar[u];
        }
    }

    x->seq.s[qlen] = '\0';
    x->seq.n = qlen;
}


static void seqenc_decode_reference_alignment(seqenc_t* E, short_read_t* r, size_t seqlen)
{
    const twobit_t* refseq = seqmap_get(E->ref, (const char*) r->seqname.s);

    if (refseq == NULL) {
        quip_error(
            "A read was aligned to sequence %s, which was not found in the reference.",
            r->seqname.s);
    }

    str_reserve(&r->seq, seqlen + 1);
    r->seq.n = 0;

    /* decode N mask */
    size_t i;
    reserve_nmask(E, seqlen);
    memset(r->seq.s, '\0', seqlen + 1);
    for (i = 0; i < seqlen; ++i) {
        if (dist2_decode(E->ac, &E->d_nmask[i])) r->seq.s[i] = 'N';
    }

    uint32_t ref_pos   = r->pos;
    uint32_t read_pos  = 0;

    i = 0;    /* cigar operation */
    size_t j; /* position within the cigar op */

    kmer_t y; /* reference nucleotide */

    for (i = 0; i < r->cigar.n; ++i) {
        switch (r->cigar.ops[i]) {
            case BAM_CEQUAL:
            case BAM_CDIFF:
            case BAM_CMATCH:
                for (j = 0; j < r->cigar.lens[i]; ++j, ++read_pos, ++ref_pos) {
                    if (r->seq.s[read_pos] == 'N') continue;

                    if (dist2_decode(E->ac, &E->d_ref_match) == SEQENC_REF_MATCH) {
                        y = twobit_get(refseq, ref_pos);
                        r->seq.s[read_pos] = kmertochar[y];
                    }
                    else {
                        r->seq.s[read_pos] = kmertochar[dist4_decode(E->ac, &E->d_ref_ins_nuc)];
                    }
                }
                break;

            case BAM_CINS:
                for (j = 0; j < r->cigar.lens[i]; ++j, ++read_pos) {
                    if (r->seq.s[read_pos] == 'N') continue;
                    r->seq.s[read_pos] = kmertochar[dist4_decode(E->ac, &E->d_ref_ins_nuc)];
                }
                break;

            case BAM_CDEL:
                ref_pos += r->cigar.lens[i];
                break;

            case BAM_CREF_SKIP:
                ref_pos += r->cigar.lens[i];
                break;

            case BAM_CSOFT_CLIP:
                for (j = 0; j < r->cigar.lens[i]; ++j, ++read_pos) {
                    if (r->seq.s[read_pos] == 'N') continue;
                    r->seq.s[read_pos] = kmertochar[dist4_decode(E->ac, &E->d_ref_ins_nuc)];
                }
                break;

            case BAM_CHARD_CLIP:
                ref_pos += r->cigar.lens[i];
                break;

            case BAM_CPAD:
                quip_error("Unsupported cigar operation.");
                break;
        }
    }
    r->seq.s[seqlen] = '\0';
    r->seq.n = seqlen;

    if (read_pos != seqlen) {
        quip_error("Cigar operations do not account for full read length.");
    }

    if (r->strand) str_revcomp(r->seq.s, r->seq.n);
}


void seqenc_decode(seqenc_t* E, short_read_t* x, size_t n)
{
    if (E->ref != NULL && (x->flags & BAM_FUNMAP) == 0) {
        seqenc_decode_reference_alignment(E, x, n);
    }
    else {
        uint32_t type = dist2_decode(E->ac, &E->d_type);

        if (type == SEQENC_TYPE_SEQUENCE) seqenc_decode_seq(E, x, n);
        else                              seqenc_decode_alignment(E, x, n);
    }
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


