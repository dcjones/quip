 
#include "assembler.h"
#include "bloom.h"
#include "kmer.h"
#include "kmerhash.h"
#include "misc.h"
#include "seqenc.h"
#include "seqset.h"
#include "twobit.h"
#include <assert.h>
#include <limits.h>
#include <string.h>
#include <math.h>


/* Maximum score for an alignment to be reported, in proportion
 * of positions that mismatch. */
static const double max_align_score = 0.22;

/* Minimum allowable contig length. */
static size_t min_contig_len = 200;


static uint32_t read_uint32(quip_reader_t reader, void* reader_data)
{
    uint8_t bytes[4];
    size_t cnt = reader(reader_data, bytes, 4);

    if (cnt < 4) {
        fprintf(stderr, "Unexpected end of file.\n");
        exit(EXIT_FAILURE);
    }

    return ((uint32_t) bytes[0] << 24) |
           ((uint32_t) bytes[1] << 16) |
           ((uint32_t) bytes[2] << 8) |
           ((uint32_t) bytes[3]);
}


static void write_uint32(quip_writer_t writer, void* writer_data, uint32_t x)
{
    uint8_t bytes[4] = { (uint8_t) (x >> 24),
                         (uint8_t) (x >> 16),
                         (uint8_t) (x >> 8),
                         (uint8_t) x };
    writer(writer_data, bytes, 4);
}


struct assembler_t_
{
    /* don't actually assemble anything */
    bool quick;

    /* function to write compressed output */
    quip_writer_t writer;
    void* writer_data;

    /* kmer bit-mask used in assembly */
    kmer_t assemble_kmer_mask;

    /* kmer bit-mask used in alignment */
    kmer_t align_kmer_mask;

    /* read set */
    seqset_t* S;

    /* reads, in-order */
    uint32_t* ord;

    /* size of the ord array */
    size_t ord_size;

    /* total number of reads */
    uint32_t N;

    /* current read */
    twobit_t* x;

    /* k-mer size used for assembly */
    size_t assemble_k;

    /* k-mer size used for seeds of alignment */
    size_t align_k;

    /* k-mer table used for assembly */
    bloom_t* B;

    /* k-mer table used for alignment */
    kmerhash_t* H;

    /* count required to be nominated as a seed candidate */
    unsigned int count_cutoff;

    /* nucleotide sequence encoder */
    seqenc_t* seqenc;

    /* assembled contigs and their reverse complements */
    twobit_t* supercontig;
    twobit_t* supercontig_rc;

    /* should we try attempt to assemble contigs with the current batch of reads
     * */
    bool assemble_batch;
};


static void build_kmer_hash(assembler_t* A);


assembler_t* assembler_alloc(
        quip_writer_t writer, void* writer_data,
        size_t assemble_k, size_t align_k, bool quick)
{
    assembler_t* A = malloc_or_die(sizeof(assembler_t));
    memset(A, 0, sizeof(assembler_t));

    A->quick = quick;

    A->writer = writer;
    A->writer_data = writer_data;

    A->assemble_k = assemble_k;
    A->align_k    = align_k;

    A->assemble_kmer_mask = 0;
    size_t i;
    for (i = 0; i < A->assemble_k; ++i) {
        A->assemble_kmer_mask = (A->assemble_kmer_mask << 2) | 0x3;
    }

    A->align_kmer_mask = 0;
    for (i = 0; i < A->align_k; ++i) {
        A->align_kmer_mask = (A->align_kmer_mask << 2) | 0x3;
    }

    /* We delay allocation of the sequence encoder. To keep memory usage under
     * control we first finish assembling and free the bloom filter before
     * allocating seqenc. */
    A->seqenc = NULL;

    A->assemble_batch = true;

    /* If we are not assembling, we do not need any of the data structure
     * initialized below. */
    if (!quick) {
        A->S = seqset_alloc();

        A->ord_size = 1024;
        A->ord = malloc_or_die(A->ord_size * sizeof(const twobit_t*));

        A->N = 0;

        A->B = bloom_alloc(4097152, 8);

        A->x = twobit_alloc();

        A->H = kmerhash_alloc();

        A->count_cutoff = 2;
    }
    else {
        A->seqenc = seqenc_alloc_encoder(writer, writer_data);
    }

    return A;
}


void assembler_free(assembler_t* A)
{
    seqset_free(A->S);
    free(A->ord);
    bloom_free(A->B);
    kmerhash_free(A->H);
    twobit_free(A->x);
    seqenc_free(A->seqenc);
    twobit_free(A->supercontig);
    twobit_free(A->supercontig_rc);
    free(A);
}




/* Try desperately to align the given read. If no good alignment is found, just
 * encode the sequence. Return 0 if no alignment was found, 1 if a seed was
 * found but no good alignment, and 2 if a good alignment was found. */
static bool align_read(assembler_t* A, const twobit_t* seq)
{
    /* We only consider the first few seed hits found in the hash table. The
     * should be in an approximately random order. */
    static const size_t max_seeds = 100;

    /* position of the seed with the subject and query sequecne, resp. */
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

    /* Don't try to align any reads that are shorer than the seed length */
    qlen = twobit_len(seq);
    if ((size_t) qlen < A->align_k) {
        seqenc_encode_twobit_seq(A->seqenc, seq);
        return false;
    }

    uint32_t max_mismatch = 1 + max_align_score * (double) qlen;
    size_t i, j;
    for (i = 0; i < 3 && best_aln_score > 0.0; ++i) {
        if      (i == 0) qpos = 0;
        else if (i == 1) qpos = 2 * A->align_k; 
        else if (i == 2) qpos = 4 * A->align_k;
        if (qpos + A->align_k > (size_t) qlen) qpos = qlen - A->align_k - 1;


        x = twobit_get_kmer_rev(
                seq,
                qpos,
                A->align_k);
        y = kmer_canonical(x, A->align_k);

        poslen = kmerhash_get(A->H, y, &pos);
        poslen = poslen > max_seeds ? max_seeds : poslen;


        for (j = 0; j < poslen && best_aln_score > 0.0; ++j, ++pos) {

            if (*pos >= 0) {
                if (x == y) {
                    spos = *pos;
                    strand = 0;
                }
                else {
                    spos = slen - *pos - A->align_k;
                    strand = 1;
                }
            }
            else {
                if (x == y) {
                    spos = slen + *pos - A->align_k + 1;
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
            A->seqenc, best_spos, best_strand, seq);

        return true;
    }
    else {
        seqenc_encode_twobit_seq(A->seqenc, seq);
        return false;
    }
}



void assembler_add_seq(assembler_t* A, const char* seq, size_t seqlen)
{
    if (A->quick) {
        seqenc_encode_char_seq(A->seqenc, seq, seqlen);
    }
    else if (A->assemble_batch) {
        if (A->N == A->ord_size) {
            A->ord_size *= 2;
            A->ord = realloc_or_die(A->ord, A->ord_size * sizeof(uint32_t));
        }

        /* does the read contain non-nucleotide characters ? */
        size_t i;
        bool has_N = false;
        for (i = 0; i < seqlen; ++i) {
            if (seq[i] == 'N') {
                has_N = true;
                break;
            }
        }

        if (has_N) {
            A->ord[A->N++] = seqset_inc_eb(A->S, seq);
        }
        else {
            twobit_copy_n(A->x, seq, seqlen);
            A->ord[A->N++] = seqset_inc_tb(A->S, A->x);
        }
    }
    else {
        /* does the read contain non-nucleotide characters ? */
        size_t i;
        bool has_N = false;
        for (i = 0; i < seqlen; ++i) {
            if (seq[i] == 'N') {
                has_N = true;
                break;
            }
        }

        if (has_N) {
            seqenc_encode_char_seq(A->seqenc, seq, seqlen);
        }
        else {
            twobit_copy_n(A->x, seq, seqlen);
            align_read(A, A->x);
        }

    }
}



static void make_contig(assembler_t* A, twobit_t* seed, twobit_t* contig)
{
    twobit_clear(contig);


    /* delete all kmers in the seed */
    kmer_t x = twobit_get_kmer_rev(seed, 0, A->assemble_k);
    size_t i;
    for (i = A->assemble_k; i < twobit_len(seed); ++i) {
        bloom_ldec(A->B, kmer_canonical((x << 2) | twobit_get(seed, i), A->assemble_k));
    }


    /* expand the contig as far left as possible */
    unsigned int cnt, cnt2, cnt_best, cnt2_best;

    kmer_t nt, nt2, nt_best = 0, xc, y, z;


    /* Greedily append nucleotides to the contig using 
     * approximate k-mer counts stored in the bloom filter. */
    x = twobit_get_kmer_rev(seed, 0, A->assemble_k);
    while (true) {
        bloom_ldec(A->B, kmer_canonical(x, A->assemble_k));

        x = (x >> 2) & A->assemble_kmer_mask;
        cnt_best = 0;
        for (nt = 0; nt < 4; ++nt) {
            y = nt << (2 * (A->assemble_k - 1));
            xc = kmer_canonical(x | y, A->assemble_k);
            cnt = bloom_get(A->B, xc);

            /* Look ahead two k-mers for a somewhat better
             * greedy choice. */
            cnt2_best = 0;
            for (nt2 = 0; nt2 < 4; ++nt2) {
                z = (nt2 << (2 * (A->assemble_k - 1))) |
                    (nt  << (2 * (A->assemble_k - 2))) |
                    (x >> 2);
                z &= A->assemble_kmer_mask; 
                cnt2 = bloom_get(A->B, kmer_canonical(z, A->assemble_k));
                if (cnt2 > cnt2_best) cnt2_best = cnt2;
            }

            if (cnt + cnt2_best > cnt_best) {
                cnt_best = cnt + cnt2_best;
                nt_best  = nt;
            }
        }

        if (cnt_best > 0) {
            y = nt_best << (2 * (A->assemble_k - 1));
            x = x | y;
            twobit_append_kmer(contig, nt_best, 1);
        }
        else break;
    }

    twobit_reverse(contig);
    twobit_append_twobit(contig, seed);

    x = twobit_get_kmer_rev(seed, twobit_len(seed) - A->assemble_k, A->assemble_k);
    while (true) {
        bloom_ldec(A->B, kmer_canonical(x, A->assemble_k));

        x = (x << 2) & A->assemble_kmer_mask;
        cnt_best = 0;
        for (nt = 0; nt < 4; ++nt) {
            xc = kmer_canonical(x | nt, A->assemble_k);
            cnt = bloom_get(A->B, xc);

            cnt2_best = 0;
            for (nt2 = 0; nt2 < 4; ++nt2) {
                z = (x << 2) | (nt << 2) | nt2;
                z &= A->assemble_kmer_mask;
                cnt2 = bloom_get(A->B, kmer_canonical(z, A->assemble_k));
                if (cnt2 > cnt2_best) cnt2_best = cnt2;
            }

            if (cnt + cnt2_best > cnt_best) {
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


static int seqset_value_cnt_cmp(const void* a_, const void* b_)
{
    seqset_value_t* a = (seqset_value_t*) a_;
    seqset_value_t* b = (seqset_value_t*) b_;

    if      (a->cnt < b->cnt) return 1;
    else if (a->cnt > b->cnt) return -1;
    /* break ties with the index, just to make things deterministic */
    else if (a->idx < b->idx) return -1;
    else if (a->idx > b->idx) return 1;
    else                      return 0;
}


static int seqset_value_idx_cmp(const void* a_, const void* b_)
{
    seqset_value_t* a = (seqset_value_t*) a_;
    seqset_value_t* b = (seqset_value_t*) b_;

    if      (a->idx < b->idx) return -1;
    else if (a->idx > b->idx) return 1;
    else                      return 0;
}


static void build_kmer_hash(assembler_t* A)
{
    kmerhash_clear(A->H);

    size_t len = twobit_len(A->supercontig);
    size_t pos;
    kmer_t x = 0, y;
    for (pos = 0; pos < len; ++pos) {
        x = ((x << 2) | twobit_get(A->supercontig, pos)) & A->align_kmer_mask;

        if (pos + 1 >= A->align_k) {
            y = kmer_canonical(x, A->align_k);
            if (x == y) kmerhash_put(A->H, y, pos + 1 - A->align_k);
            else        kmerhash_put(A->H, y, - (int32_t) (pos + 2 - A->align_k));
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


/* Align buffered reads to contigs. */
static void align_to_contigs(assembler_t* A,
                             seqset_value_t* xs, size_t xs_len)
{
    if (quip_verbose) fprintf(stderr, "aligning reads to contigs ...\n");

    /* Sort reads by index. Afterwards xs[i].idx == i should be true for every
     * i, where 0 <= i < xs_len */
    qsort(xs, xs_len, sizeof(seqset_value_t), seqset_value_idx_cmp);
    const char* ebseq;
    size_t aln_count = 0;
    size_t aborted_aln_count = 0;
    size_t i;
    int ret;

    for (i = 0; i < A->N; ++i) {
        if (xs[A->ord[i]].is_twobit) {
            ret = align_read(A, xs[A->ord[i]].seq.tb);
            if (ret == 1)      ++aln_count;
        }
        else {
            ebseq = xs[A->ord[i]].seq.eb;
            seqenc_encode_char_seq(A->seqenc, ebseq, strlen(ebseq));
        }
    }

    if (quip_verbose) {
        fprintf(stderr,
                "\t%zu / %zu [%0.2f%%] reads aligned to contigs (%0.2f%% aborted)\n",
                aln_count, (size_t) A->N,
                100.0 * (double) aln_count / (double) A->N,
                100.8 * (double) aborted_aln_count / (double) A->N);

        fprintf(stderr, "done.\n");
    }

}



/* Count the number of occurrences of each k-mer in the array of reads xs, of
 * length n. */
static void count_kmers(assembler_t* A, seqset_value_t* xs, size_t n)
{
    if (quip_verbose) fprintf(stderr, "counting k-mers ... ");

    size_t i, j, seqlen;
    kmer_t x, y;

    for (i = 0; i < n; ++i) {
        if (!xs[i].is_twobit) continue;
        seqlen = twobit_len(xs[i].seq.tb);
        x = 0;
        for (j = 0; j < seqlen; ++j) {
            x = ((x << 2) | twobit_get(xs[i].seq.tb, j)) & A->assemble_kmer_mask;

            if (j + 1 >= A->assemble_k) {
                y = kmer_canonical(x, A->assemble_k);
                bloom_add(A->B, y, xs[i].cnt);
            }
        }
    }

    if (quip_verbose) fprintf(stderr, "done.\n");
}


/* Heuristically build contigs from k-mers counts and a set of reads */
static void make_contigs(assembler_t* A, seqset_value_t* xs, size_t n)
{
    if (quip_verbose) fprintf(stderr, "assembling contigs ... ");

    A->supercontig = twobit_alloc();
    twobit_t* contig = twobit_alloc();
    size_t len;

    size_t contig_cnt = 0;

    size_t i;
    for (i = 0; i < n && xs[i].cnt >= A->count_cutoff; ++i) {
        if (!xs[i].is_twobit) continue;

        make_contig(A, xs[i].seq.tb, contig);

        /* reject over terribly short contigs */
        len = twobit_len(contig);
        if (len < min_contig_len) {
            continue;
        }

        twobit_append_twobit(A->supercontig, contig);
        ++contig_cnt;
    }

    twobit_free(contig);

    A->supercontig_rc = twobit_alloc_n(twobit_len(A->supercontig));
    twobit_revcomp(A->supercontig_rc, A->supercontig);

    if (quip_verbose) fprintf(stderr, "done. (%zu contigs, %zunt)\n",
                              contig_cnt, twobit_len(A->supercontig));
}



size_t assembler_finish(assembler_t* A)
{
    size_t bytes = 0;

    if (!A->quick && A->assemble_batch) {

        /* dump reads and sort by copy number*/
        seqset_value_t* xs = seqset_dump(A->S);
        qsort(xs, seqset_size(A->S), sizeof(seqset_value_t), seqset_value_cnt_cmp);

        size_t n = seqset_size(A->S);

        count_kmers(A, xs, n);
        make_contigs(A, xs, n);

        /* Free up memory we won't need any more */
        bloom_free(A->B);
        A->B = NULL;

        /* Now that we have some memory to spare, bring up the sequence encoder. */
        A->seqenc = seqenc_alloc_encoder(A->writer, A->writer_data);
        seqenc_set_supercontig(A->seqenc, A->supercontig);

        /* Number of bytes needed to write the contig count and contig lengths */
        bytes += 4;

        /* write the supercontig */
        seqenc_encode_twobit_seq(A->seqenc, A->supercontig);

        index_contigs(A);

        align_to_contigs(A, xs, n);

        free(xs);
        seqset_free(A->S);
        A->S = NULL;

        free(A->ord);
        A->ord = NULL;
    }
    else {
        /* Zero contigs in the block. */
        bytes += 4;
    }

    bytes += seqenc_finish(A->seqenc);

    if (!A->quick) {
        seqenc_get_supercontig_consensus(A->seqenc, A->supercontig);
        index_contigs(A);
    }

    return bytes;
}


void assembler_flush(assembler_t* A)
{
    /* write supercontig length */
    if (!A->quick && A->assemble_batch) {
        write_uint32(A->writer, A->writer_data, twobit_len(A->supercontig));

        /* only assemble the first block. */
        A->assemble_batch = false;
    }
    else {
        write_uint32(A->writer, A->writer_data, 0);
    }

    seqenc_flush(A->seqenc);
}



struct disassembler_t_
{
    seqenc_t* seqenc;

    quip_reader_t reader;
    void* reader_data;

    bool init_state;
};


disassembler_t* disassembler_alloc(quip_reader_t reader, void* reader_data)
{
    disassembler_t* D = malloc_or_die(sizeof(disassembler_t));

    D->seqenc = seqenc_alloc_decoder(reader, reader_data);
    D->reader = reader;
    D->reader_data = reader_data;
    D->init_state = true;

    return D;
}


void disassembler_free(disassembler_t* D)
{
    if (D == NULL) return;
    seqenc_free(D->seqenc);
    free(D);
}


void disassembler_read(disassembler_t* D, seq_t* x, size_t n)
{
    if (D->init_state) {
        uint32_t supercontig_len = read_uint32(D->reader, D->reader_data);

        if (supercontig_len > 0) {
            seqenc_prepare_decoder(D->seqenc, supercontig_len);
        }

        seqenc_start_decoder(D->seqenc);

        D->init_state = false;
    }

    seqenc_decode(D->seqenc, x, n);
}


void disassembler_reset(disassembler_t* D)
{
    if (!D->init_state) seqenc_reset_decoder(D->seqenc);
    D->init_state = true;
}



