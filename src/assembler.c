
#include "assembler.h"
#include "bloom.h"
#include "kmer.h"
#include "kmerhash.h"
#include "misc.h"
#include "seqenc.h"
#include "twobit.h"
#include "sam/bam.h"
#include <assert.h>
#include <limits.h>
#include <string.h>
#include <math.h>
#include <inttypes.h>

/* k-mer used for de bruijn graph assembly */
static const size_t assemble_k = 25;

/* bitmask for 2-bit encoded k-mers */
static const kmer_t assemble_kmer_mask = 0x0003fffffffffffful;

/* k-mer size used for seed-and-extend alignment */
static const size_t align_k = 12;

/* bitmask for 2-bit encoded k-mers */
static const kmer_t align_kmer_mask = 0x0000000000fffffful;

/* Maximum score for an alignment to be reported, in proportion
 * of positions that mismatch. */
static const double max_align_score = 0.22;

/* Minimum allowable contig length. */
static const size_t min_contig_len = 200;

/* Sizes of fixed-size data structures. */
static const size_t seeds_n = 200000;
static const size_t bloom_n = 4000000;
static const size_t bloom_m = 8;

/* Number of reads to be used for assembly. */
size_t quip_assembly_n = 2500000;


struct assembler_t_
{
    /* actually assemble something */
    bool assemble;

    /* quip version field */
    uint8_t quip_version;

    /* function to write compressed output */
    quip_writer_t writer;
    void* writer_data;

    /* candidate seeds */
    twobit_t** seeds;
    size_t seeds_len;

    /* current read */
    twobit_t* x;

    /* k-mer table used for assembly */
    bloom_t* B;

    /* k-mer table used for alignment */
    kmerhash_t* H;

    /* nucleotide sequence encoder */
    seqenc_t* seqenc;

    /* assembled contigs and their reverse complements */
    twobit_t* supercontig;
    twobit_t* supercontig_rc;

    /* reference, for reference based alignment */
    const seqmap_t* ref;

    /* Number of reads before assembly is triggered. */
    size_t assembly_pending_n;

    /* statistics used in verbose reporting */
    uint64_t stat_n;
    uint64_t stat_aligned_count;
    uint64_t stat_assemble_count;
};


static void build_kmer_hash(assembler_t* A);


assembler_t* assembler_alloc(
        quip_writer_t   writer,
        void*           writer_data,
        bool            assemble,
        uint8_t         quip_version,
        const seqmap_t* ref)
{
    assembler_t* A = malloc_or_die(sizeof(assembler_t));
    memset(A, 0, sizeof(assembler_t));

    A->assemble = assemble;
    A->quip_version = quip_version;

    A->writer = writer;
    A->writer_data = writer_data;
    A->ref = ref;

    /* If we are not assembling, we do not need any of the data structure
     * initialized below. */
    if (assemble) {
        A->seeds = malloc_or_die(seeds_n * sizeof(twobit_t*));
        memset(A->seeds, 0, seeds_n * sizeof(twobit_t*));
        A->seeds_len = 0;

        A->B = bloom_alloc(bloom_n, bloom_m);
        A->x = twobit_alloc();
        A->H = kmerhash_alloc();
        A->assembly_pending_n = quip_assembly_n;
    }

    A->seqenc = seqenc_alloc_encoder(writer, writer_data, ref);

    return A;
}


void assembler_free(assembler_t* A)
{
    if (A == NULL) return;

    size_t i;
    for (i = 0; i < A->seeds_len; ++i) {
        twobit_free(A->seeds[i]);
    }
    free(A->seeds);

    bloom_free(A->B);
    kmerhash_free(A->H);
    twobit_free(A->x);
    seqenc_free(A->seqenc);
    twobit_free(A->supercontig);
    twobit_free(A->supercontig_rc);
    free(A);
}


/* Count the number of occurrences of each k-mer in a read x */
static void count_kmers(bloom_t* B, const twobit_t* seq)
{
    size_t i, seqlen;
    kmer_t x, y;

    seqlen = twobit_len(seq);
    x = 0;
    for (i = 0; i < seqlen; ++i) {
        x = ((x << 2) | twobit_get(seq, i)) & assemble_kmer_mask;

        if (i + 1 >= assemble_k) {
            y = kmer_canonical(x, assemble_k);
            bloom_add(B, y, 1);
        }
    }
}


/* Try desperately to align the given read. If no good alignment is found, just
 * encode the sequence. Return true if an alignment was found. */
static bool align_read(assembler_t* A, const unsigned char* seq_str, const twobit_t* seq)
{
    /* We only consider the first few seed hits found in the hash table. The
     * should be in an approximately random order. */
    static const size_t max_seeds = 100;

    /* position of the seed with the subject and query sequence, resp. */
    int spos, qpos;

    /* subject and query lengths, reps */
    int qlen, slen = twobit_len(A->supercontig);
    kmer_t x, y;
    uint8_t strand;

    /* positions matching the seed k-mer */
    kmer_pos_t* pos;
    size_t poslen;

    twobit_t* contig;

    /* optimal alignment found so far */
    double   best_aln_score = HUGE_VAL;
    uint32_t best_spos = 0;
    uint8_t  best_strand = 0;

    double aln_score;

    /* Don't try to align any reads that are shorter than the seed length */
    qlen = twobit_len(seq);
    if ((size_t) qlen < align_k) {
        seqenc_encode_twobit_seq(A->seqenc, seq_str, seq);
        return false;
    }

    uint32_t max_mismatch = 1 + max_align_score * (double) qlen;
    size_t i, j;
    for (i = 0; i < 3 && best_aln_score > 0.0; ++i) {
        if      (i == 0) qpos = 0;
        else if (i == 1) qpos = 2 * align_k;
        else if (i == 2) qpos = 4 * align_k;
        if (qpos + align_k > (size_t) qlen) qpos = qlen - align_k - 1;


        x = twobit_get_kmer_rev(
                seq,
                qpos,
                align_k);
        y = kmer_canonical(x, align_k);

        poslen = kmerhash_get(A->H, y, &pos);
        poslen = poslen > max_seeds ? max_seeds : poslen;


        for (j = 0; j < poslen && best_aln_score > 0.0; ++j, ++pos) {

            if (*pos >= 0) {
                if (x == y) {
                    spos = *pos;
                    strand = 0;
                }
                else {
                    spos = slen - *pos - align_k;
                    strand = 1;
                }
            }
            else {
                if (x == y) {
                    spos = slen + *pos - align_k + 1;
                    strand = 1;
                }
                else {
                    spos = -*pos - 1;
                    strand = 0;
                }
            }

            /* Is a full alignment possible with this seed? */
            if (spos < qpos || slen - spos < qlen - qpos) {
                continue;
            }

            contig = strand == 0 ? A->supercontig : A->supercontig_rc;

            aln_score = (double) twobit_mismatch_count(contig, seq, spos - qpos, max_mismatch) / (double) qlen;

            if (aln_score <= max_align_score &&
                aln_score < best_aln_score)
            {
                best_aln_score  = aln_score;
                best_strand     = strand;
                best_spos       = spos - qpos;
            }
        }
    }

    if (best_aln_score < HUGE_VAL) {
        seqenc_encode_alignment(
            A->seqenc, best_spos, best_strand, seq_str, seq);

        return true;
    }
    else {
        seqenc_encode_twobit_seq(A->seqenc, seq_str, seq);
        return false;
    }
}


static void make_contig(uint8_t quip_version, bloom_t* B,
                        twobit_t* seed, twobit_t* contig)
{
    twobit_clear(contig);

    /* delete all kmers in the seed */
    kmer_t x = twobit_get_kmer_rev(seed, 0, assemble_k);
    size_t i;
    for (i = assemble_k; i < twobit_len(seed); ++i) {
        bloom_ldec(B, kmer_canonical((x << 2) | twobit_get(seed, i), assemble_k));
    }

    /* expand the contig as far left as possible */
    unsigned int cnt, cnt2, cnt_best, cnt2_best;

    kmer_t nt, nt2, nt_best = 0, xc, y, z;

    /* Greedily append nucleotides to the contig using
     * approximate k-mer counts stored in the bloom filter. */
    x = twobit_get_kmer_rev(seed, 0, assemble_k);
    while (true) {
        bloom_ldec(B, kmer_canonical(x, assemble_k));

        x = (x >> 2) & assemble_kmer_mask;
        cnt_best = 0;
        for (nt = 0; nt < 4; ++nt) {
            y = nt << (2 * (assemble_k - 1));
            xc = kmer_canonical(x | y, assemble_k);
            cnt = bloom_get(B, xc);

            /* Look ahead two k-mers for a somewhat better
             * greedy choice. */
            cnt2_best = 0;
            for (nt2 = 0; nt2 < 4; ++nt2) {
                z = (nt2 << (2 * (assemble_k - 1))) |
                    (nt  << (2 * (assemble_k - 2))) |
                    (x >> 2);
                z &= assemble_kmer_mask;
                cnt2 = bloom_get(B, kmer_canonical(z, assemble_k));
                if (cnt2 > cnt2_best) cnt2_best = cnt2;
            }

            if ((quip_version == 2 && cnt + cnt2_best > cnt_best) ||
                (quip_version == 3 && cnt > 0 && cnt + cnt2_best > cnt_best)) {
                cnt_best = cnt + cnt2_best;
                nt_best  = nt;
            }
        }

        if (cnt_best > 0) {
            y = nt_best << (2 * (assemble_k - 1));
            x = x | y;
            twobit_append_kmer(contig, nt_best, 1);
        }
        else break;
    }

    twobit_reverse(contig);
    twobit_append_twobit(contig, seed);

    x = twobit_get_kmer_rev(seed, twobit_len(seed) - assemble_k, assemble_k);
    while (true) {
        bloom_ldec(B, kmer_canonical(x, assemble_k));

        x = (x << 2) & assemble_kmer_mask;
        cnt_best = 0;
        for (nt = 0; nt < 4; ++nt) {
            xc = kmer_canonical(x | nt, assemble_k);
            cnt = bloom_get(B, xc);

            cnt2_best = 0;
            for (nt2 = 0; nt2 < 4; ++nt2) {
                z = (x << 2) | (nt << 2) | nt2;
                z &= assemble_kmer_mask;
                cnt2 = bloom_get(B, kmer_canonical(z, assemble_k));
                if (cnt2 > cnt2_best) cnt2_best = cnt2;
            }

            if ((quip_version == 2 && cnt + cnt2_best > cnt_best) ||
                (quip_version == 3 && cnt > 0 && cnt + cnt2_best > cnt_best)) {
                cnt_best = cnt + cnt2_best;
                nt_best  = nt;
            }
        }

        if (cnt_best > 0) {
            x = x | nt_best;
            twobit_append_kmer(contig, nt_best, 1);
        }
        else break;
    }
}



/* Heuristically build contigs from k-mers counts and a set of reads */
static void make_contigs(
    uint8_t quip_version,
    twobit_t** supercontig,
    twobit_t** supercontig_rc,
    bloom_t* B,
    twobit_t** seeds, size_t n)
{
    if (quip_verbose) fprintf(stderr, "assembling contigs ...\n");

    *supercontig = twobit_alloc();
    twobit_t* contig = twobit_alloc();
    size_t len;
    size_t i;

    size_t contig_cnt = 0;

    for (i = 0; i < n; ++i) {
        make_contig(quip_version, B, seeds[i], contig);

        /* reject over terribly short contigs */
        len = twobit_len(contig);
        if (len < min_contig_len) {
            continue;
        }

        twobit_append_twobit(*supercontig, contig);
        ++contig_cnt;
    }

    twobit_free(contig);

    *supercontig_rc = twobit_alloc_n(twobit_len(*supercontig));
    twobit_revcomp(*supercontig_rc, *supercontig);

    if (quip_verbose) fprintf(stderr, "done. (%zu contigs, %zunt)\n",
                              contig_cnt, twobit_len(*supercontig));
}

static void build_kmer_hash(assembler_t* A)
{
    kmerhash_clear(A->H);

    size_t len = twobit_len(A->supercontig);
    size_t pos;
    kmer_t x = 0, y;
    for (pos = 0; pos < len; ++pos) {
        x = ((x << 2) | twobit_get(A->supercontig, pos)) & align_kmer_mask;

        if (pos + 1 >= align_k) {
            y = kmer_canonical(x, align_k);
            if (x == y) kmerhash_put(A->H, y, pos + 1 - align_k);
            else        kmerhash_put(A->H, y, - (int32_t) (pos + 2 - align_k));
        }
    }
}


/* build new indexes after updating the consensus sequences */
static void index_contigs(assembler_t* A)
{
    if (quip_verbose) fprintf(stderr, "indexing contigs ... ");

    build_kmer_hash(A);
    twobit_revcomp(A->supercontig_rc, A->supercontig);

    if (quip_verbose) fprintf(stderr, "done.\n");
}


void assembler_add_seq(assembler_t* A, const short_read_t* seq)
{
    A->stat_n++;

    seqenc_encode_extras(A->seqenc, seq);

    if (A->ref != NULL && (seq->flags & BAM_FUNMAP) == 0) {
        seqenc_encode_reference_alignment(A->seqenc, seq);
        A->stat_aligned_count++;
        return;
    }

    if (!A->assemble) {
        seqenc_encode_char_seq(A->seqenc, seq->seq.s, seq->seq.n);
    }
    else if (A->assembly_pending_n > 0) {
        twobit_copy_str_n(A->x, (char*) seq->seq.s, seq->seq.n);

        if (A->assembly_pending_n <= seeds_n - A->seeds_len) {
            A->seeds[A->seeds_len] =
                twobit_alloc_n(seq->seq.n);
            twobit_copy(A->seeds[A->seeds_len], A->x);
            A->seeds_len++;
        }

        count_kmers(A->B, A->x);
        seqenc_encode_twobit_seq(A->seqenc, seq->seq.s, A->x);
        --A->assembly_pending_n;

        if (A->assembly_pending_n == 0) {
            make_contigs(A->quip_version, &A->supercontig, &A->supercontig_rc,
                         A->B, A->seeds, A->seeds_len);

            /* Free up memory we won't need any more */
            bloom_free(A->B);
            A->B = NULL;

            size_t i;
            for (i = 0; i < A->seeds_len; ++i) {
                twobit_free(A->seeds[i]);
            }
            free(A->seeds);
            A->seeds = 0;
            A->seeds_len = 0;

            seqenc_set_supercontig(A->seqenc, A->supercontig);

            index_contigs(A);
        }
    }
    else {
        twobit_copy_str_n(A->x, (char*) seq->seq.s, seq->seq.n);
        if (align_read(A, seq->seq.s, A->x)) A->stat_assemble_count++;
    }
}


size_t assembler_finish(assembler_t* A)
{
    size_t bytes = seqenc_finish(A->seqenc);

    /* This leads to slightly better compression, but slows things
       down a lot on large assemblies. */
/*
    if (A->assemble && A->assembly_pending_n == 0) {
        seqenc_get_supercontig_consensus(A->seqenc, A->supercontig);
        index_contigs(A);
    }
*/

    if (quip_verbose) {
        fprintf(stderr, "%2.1f%% aligned to reference.\n",
            100.0 * (double) A->stat_aligned_count / (double) A->stat_n);
        fprintf(stderr, "%2.1f%% aligned to assembled contigs.\n",
            100.0 * (double) A->stat_assemble_count / (double) A->stat_n);
    }

    A->stat_n = 0;
    A->stat_aligned_count = 0;
    A->stat_assemble_count = 0;

    return bytes;
}


void assembler_flush(assembler_t* A)
{
    seqenc_flush(A->seqenc);
}


struct disassembler_t_
{
    /* actually assemble something */
    bool assemble;

    /* quip header version used during compression */
    uint8_t quip_version;

    /* function to read compressed input */
    quip_reader_t reader;
    void* reader_data;

    /* nucleotide sequence encoder */
    seqenc_t* seqenc;

    /* assembled contigs and their reverse complements */
    twobit_t* supercontig;
    twobit_t* supercontig_rc;

    /* candidate seeds */
    twobit_t** seeds;
    size_t seeds_len;

    /* current read */
    twobit_t* x;

    /* k-mer table used for assembly */
    bloom_t* B;

    /* reference, for reference based alignment */
    const seqmap_t* ref;

    /* Number of reads before assembly is triggered. */
    size_t assembly_pending_n;

    /* Initial state. */
    bool initial_state;
};


disassembler_t* disassembler_alloc(
    quip_reader_t reader,
    void* reader_data,
    bool assemble,
    uint8_t quip_version,
    const seqmap_t* ref)
{
    disassembler_t* D = malloc_or_die(sizeof(disassembler_t));
    memset(D, 0, sizeof(disassembler_t));

    D->seqenc = seqenc_alloc_decoder(reader, reader_data, ref);
    D->reader = reader;
    D->reader_data = reader_data;
    D->ref = ref;
    D->assemble = assemble;
    D->quip_version = quip_version;
    D->initial_state = true;

    if (assemble) {
        D->seeds = malloc_or_die(seeds_n * sizeof(twobit_t*));
        memset(D->seeds, 0, seeds_n * sizeof(twobit_t*));
        D->seeds_len = 0;

        D->B = bloom_alloc(bloom_n, bloom_m);
        D->x = twobit_alloc();
        D->assembly_pending_n = quip_assembly_n;
    }

    return D;
}


void disassembler_free(disassembler_t* D)
{
    if (D == NULL) return;

    size_t i;
    for (i = 0; i < D->seeds_len; ++i) {
        twobit_free(D->seeds[i]);
    }
    free(D->seeds);

    seqenc_free(D->seqenc);
    bloom_free(D->B);
    twobit_free(D->x);
    twobit_free(D->supercontig);
    twobit_free(D->supercontig_rc);
    free(D);
}

void disassembler_read(disassembler_t* D, short_read_t* seq, size_t n)
{
    if (D->initial_state) {
        seqenc_start_decoder(D->seqenc);
        D->initial_state = false;
    }

    seqenc_decode_extras(D->seqenc, seq, n);
    seqenc_decode(D->seqenc, seq, n);

    if (D->assembly_pending_n > 0 &&
        (D->ref == NULL || seq->flags & BAM_FUNMAP) != 0) {
        twobit_copy_str_n(D->x, (char*) seq->seq.s, seq->seq.n);

        if (D->assembly_pending_n <= seeds_n - D->seeds_len) {
            D->seeds[D->seeds_len] =
                twobit_alloc_n(seq->seq.n);
            twobit_copy(D->seeds[D->seeds_len], D->x);
            D->seeds_len++;
        }

        count_kmers(D->B, D->x);
        --D->assembly_pending_n;

        if (D->assembly_pending_n == 0) {
            make_contigs(D->quip_version, &D->supercontig, &D->supercontig_rc,
                         D->B, D->seeds, D->seeds_len);

            /* Free up memory we won't need any more */
            bloom_free(D->B);
            D->B = NULL;

            size_t i;
            for (i = 0; i < D->seeds_len; ++i) {
                twobit_free(D->seeds[i]);
            }
            free(D->seeds);
            D->seeds = 0;
            D->seeds_len = 0;

            seqenc_set_supercontig(D->seqenc, D->supercontig);
        }
    }
}

void disassembler_reset(disassembler_t* D)
{
    seqenc_reset_decoder(D->seqenc);
    D->initial_state = true;
}

