
#include "assembler.h"
#include "bloom.h"
#include "kmer.h"
#include "kmerhash.h"
#include "misc.h"
#include "seqenc.h"
#include "seqset.h"
#include "sw.h"
#include "twobit.h"
#include <assert.h>
#include <limits.h>
#include <string.h>

struct assembler_t_
{
    /* function to write compressed output */
    quip_block_writer_t writer;
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
};


assembler_t* assembler_alloc(
        quip_block_writer_t writer, void* writer_data,
        size_t assemble_k, size_t align_k)
{
    assembler_t* A = malloc_or_die(sizeof(assembler_t));
    assert(A != NULL);

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

    A->S = seqset_alloc();

    A->ord_size = 1024;
    A->ord = malloc_or_die(A->ord_size * sizeof(const twobit_t*));

    A->N = 0;


    // TODO: set n and m in some principled way
    A->B = bloom_alloc(8388608, 8);

    A->x = twobit_alloc();


    A->H = kmerhash_alloc();

    A->count_cutoff = 2;

    A->seqenc = seqenc_alloc(5, writer, writer_data);

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
}



void assembler_add_seq(assembler_t* A, const char* seq, size_t seqlen)
{
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


typedef struct seed_t_
{
    kmer_t x;
    unsigned int count;
} seed_t;



static void assembler_make_contig(assembler_t* A, twobit_t* seed, twobit_t* contig)
{
    twobit_clear(contig);


    /* delete all kmers in the seed */
    kmer_t x = twobit_get_kmer(seed, 0, A->assemble_k);
    size_t i;
    for (i = A->assemble_k; i < twobit_len(seed); ++i) {
        bloom_del(A->B, kmer_canonical((x << 2) | twobit_get(seed, i), A->assemble_k));
    }


    /* expand the contig as far left as possible */
    unsigned int cnt, cnt_best = 0;

    kmer_t nt, nt_best = 0, xc, y;


    x = twobit_get_kmer(seed, 0, A->assemble_k);
    while (true) {
        bloom_del(A->B, kmer_canonical(x, A->assemble_k));

        x = (x >> 2) & A->assemble_kmer_mask;
        cnt_best = 0;
        for (nt = 0; nt < 4; ++nt) {
            y = nt << (2 * (A->assemble_k - 1));
            xc = kmer_canonical(x | y, A->assemble_k);
            cnt = bloom_get(A->B, xc);

            if (cnt > cnt_best) {
                cnt_best = cnt;
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

    x = twobit_get_kmer(seed, twobit_len(seed) - A->assemble_k, A->assemble_k);
    while (true) {
        bloom_del(A->B, kmer_canonical(x, A->assemble_k));

        x = (x << 2) & A->assemble_kmer_mask;
        cnt_best = 0;
        for (nt = 0; nt < 4; ++nt) {
            xc = kmer_canonical(x | nt, A->assemble_k);
            cnt = bloom_get(A->B, xc);

            if (cnt > cnt_best) {
                cnt_best = cnt;
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


static void index_contigs(assembler_t* A, twobit_t** contigs, size_t n)
{
    fprintf(stderr, "indexing contigs ... ");
    size_t i;
    size_t len;
    size_t pos;
    kmer_t x, y;
    twobit_t* contig;
    for (i = 0; i < n; ++i) {
        contig = contigs[i];
        len = twobit_len(contig);
        x = 0;
        for (pos = 0; pos < len; ++pos) {

            x = ((x << 2) | twobit_get(contig, pos)) & A->align_kmer_mask;

            if (pos + 1 >= A->align_k) {
                y = kmer_canonical(x, A->align_k);
                if (x == y) kmerhash_put(A->H, y, i, pos + 1 - A->align_k);
                else        kmerhash_put(A->H, y, i, - (int32_t) (pos + 2 - A->align_k));
            }
        }
    }

    fprintf(stderr, "done.\n");
}


static void align_to_contigs(assembler_t* A,
                             twobit_t** contigs, size_t contigs_len,
                             seqset_value_t* xs, size_t xs_len)
{
    if (verbose) fprintf(stderr, "aligning reads to contigs ...\n");

    /* Sort reads by index. Afterwards xs[i].idx == i should be true for every
     * i, where 0 <= i < xs_len */
    qsort(xs, xs_len, sizeof(seqset_value_t), seqset_value_idx_cmp);

    size_t i;
    for (i = 0; i < xs_len; ++i) assert(xs[i].idx == i);


    /* First pass: hash the initial k-mer of every read agains the contigs, and
     * store all candidate alignments. */

    if (verbose) fprintf(stderr, "\thashing seeds ... ");

    typedef struct cand_t_
    {
        uint32_t seq_idx;
        int      spos;
        uint8_t  strand;
        struct cand_t_* next;
    } cand_t;

    cand_t* cand;
    cand_t** cands = malloc_or_die(contigs_len * sizeof(cand_t*));
    memset(cands, 0, contigs_len * sizeof(cand_t*));

    /* We only consider the first few position a k-mer hashes to */
    const size_t max_kmer_pos = 5;

    size_t seqlen;

    kmer_pos_t* pos;
    size_t poslen;
    kmer_t x, y;
    for (i = 0; i < xs_len; ++i) {
        if (!xs[i].is_twobit) continue;
        seqlen = twobit_len(xs[i].seq.tb);
        if (seqlen < A->align_k) continue;

        x = twobit_get_kmer(
                xs[i].seq.tb,
                /* We choose the seed from the middle of the read as that
                 * minimizes the cost of local alignment. */
                (seqlen - A->align_k) / 2,
                A->align_k);
        y = kmer_canonical(x, A->align_k);

        poslen = kmerhash_get(A->H, y, &pos);
        poslen = poslen > max_kmer_pos ? max_kmer_pos : poslen;

        while (poslen--) {
            cand = malloc_or_die(sizeof(cand_t));
            cand->seq_idx = i;

            if (pos->contig_pos >= 0) {
                if (x == y) {
                    cand->spos = pos->contig_pos;
                    cand->strand = 0;
                }
                else {
                    cand->spos = (int) twobit_len(contigs[pos->contig_idx]) - pos->contig_pos - A->align_k;
                    cand->strand = 1;
                }
            }
            else {
                if (x == y) {
                    cand->spos = (int) twobit_len(contigs[pos->contig_idx]) + pos->contig_pos - (A->align_k - 1);
                    cand->strand = 1;
                }
                else {
                    cand->spos = -pos->contig_pos - 1;
                    cand->strand = 0;
                }
            }

            cand->next = cands[pos->contig_idx];
            cands[pos->contig_idx] = cand;
        }
    }

    if (verbose) fprintf(stderr, "done.\n");



    /* Second pass: one contig at a time, perform local alignment. */


    /* TODO: prune contigs to which no reads possibly align */


    if (verbose) fprintf(stderr, "\tlocal alignment ... ");

    typedef struct align_t_
    {
        uint32_t contig_idx;
        uint8_t  strand;
        int      aln_score;
        sw_alignment_t a;
    } align_t;

    align_t* alns = malloc_or_die(xs_len * sizeof(align_t));
    for (i = 0; i < xs_len; ++i) {
        alns[i].aln_score = INT_MAX;
        memset(&alns[i].a, 0, sizeof(sw_alignment_t));
    }

    uint32_t* contig_reindex = malloc_or_die(contigs_len * sizeof(uint32_t));
    size_t j;

    sw_t* sw;
    sw_t* sw_rc;
    twobit_t* contig_rc = twobit_alloc();

    int aln_score;

    for (i = 0; i < contigs_len; ++i) {

        if (!cands[i]) continue;

        seqenc_encode_twobit_seq(A->seqenc, contigs[i]);
        contig_reindex[i] = j++;
        
        sw = sw_alloc(contigs[i]);

        twobit_revcomp(contig_rc, contigs[i]);
        sw_rc = sw_alloc(contig_rc);

        while (cands[i]) {
            cand = cands[i];
            seqlen = twobit_len(xs[cand->seq_idx].seq.tb);

            /* local alignment */
            aln_score = sw_seeded_align(
                            cand->strand == 0 ? sw : sw_rc,
                            xs[cand->seq_idx].seq.tb,
                            cand->spos,
                            (seqlen - A->align_k) / 2,
                            A->align_k);


            /* better than the current alignment ? */
            /* Note: 'aln_score < seqlen / 2' is the simplistic heuristic I am
             * using the decide when outputing the alignment would be cheaper
             * than outputing the sequence.
             */
            if (aln_score >= 0 &&
                aln_score < (int) seqlen + (int) seqlen &&
                aln_score < alns[cand->seq_idx].aln_score)
            {
                alns[cand->seq_idx].contig_idx = i;
                alns[cand->seq_idx].strand = cand->strand;
                alns[cand->seq_idx].aln_score = aln_score;
                sw_trace(cand->strand == 0 ? sw : sw_rc, &alns[cand->seq_idx].a);
            }

            cands[i] = cand->next;
            free(cand);
        }

        sw_free(sw);
        sw_free(sw_rc);
    }


    if (verbose) fprintf(stderr, "done.\n");


    /* Lastly: output compressed reads */
    size_t aln_count = 0;
    for (i = 0; i < A->N; ++i) {
        if (alns[A->ord[i]].aln_score < INT_MAX) {
            seqenc_encode_alignment(A->seqenc,
                    alns[A->ord[i]].contig_idx, alns[A->ord[i]].strand,
                    &alns[A->ord[i]].a, xs[A->ord[i]].seq.tb);
            /*seqenc_encode_twobit_seq(A->seqenc, xs[A->ord[i]].seq.tb);*/
            ++aln_count;
        }
        else {
            if (xs[A->ord[i]].is_twobit) {
                seqenc_encode_twobit_seq(A->seqenc, xs[A->ord[i]].seq.tb);
            }
            else {
                seqenc_encode_char_seq(A->seqenc, xs[A->ord[i]].seq.eb);
            }
        }
    }

    if (verbose) {
        fprintf(stderr,
                "\t%zu / %zu [%0.2f%%] reads aligned to contigs\n",
                aln_count, (size_t) A->N, 100.0 * (double) aln_count / (double) A->N);
    }


    for (i = 0; i < xs_len; ++i) {
        free(alns[i].a.ops);
    }
    free(alns);
    free(contig_reindex);

    twobit_free(contig_rc);
    if (verbose) fprintf(stderr, "done.\n");
}




void assembler_assemble(assembler_t* A)
{
    /* dump reads and sort by abundance */
    seqset_value_t* xs = seqset_dump(A->S);
    qsort(xs, seqset_size(A->S), sizeof(seqset_value_t), seqset_value_cnt_cmp);

    /* count kmers */

    if (verbose) fprintf(stderr, "counting k-mers ... ");

    size_t i, j, n, seqlen;
    kmer_t x, y;

    n = seqset_size(A->S);
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

    if (verbose) fprintf(stderr, "done.\n");


    /* assemble contigs */

    if (verbose) fprintf(stderr, "assembling contigs ... ");

    size_t contigs_size = 512;
    size_t contigs_len  = 0;
    twobit_t** contigs = malloc_or_die(contigs_size * sizeof(twobit_t*));

    twobit_t* contig = twobit_alloc();

    for (i = 0; i < n && xs[i].cnt >= A->count_cutoff; ++i) {
        if (!xs[i].is_twobit) continue;

        assembler_make_contig(A, xs[i].seq.tb, contig);
        if (twobit_len(contig) < 3 * A->assemble_k) continue;

        /* TODO: when we discard a contig, it would be nice if we could return
         * its k-mers to the bloom filter */

        if (contigs_len == contigs_size) {
            contigs_size *= 2;
            contigs = realloc_or_die(contigs, contigs_size * sizeof(twobit_t*));
        }

        contigs[contigs_len++] = twobit_dup(contig);
    }

    twobit_free(contig);

    if (verbose) fprintf(stderr, "done. (%zu contigs)\n", contigs_len);

    index_contigs(A, contigs, contigs_len);
    align_to_contigs(A, contigs, contigs_len, xs, n);

    for (i = 0; i < contigs_len; ++i) twobit_free(contigs[i]);
    free(contigs);

    seqenc_flush(A->seqenc);
}


