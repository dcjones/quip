
#include "assembler.h"
#include "bloom.h"
#include "kmer.h"
#include "kmerhash.h"
#include "misc.h"
#include "seqset.h"
#include "sw.h"
#include "twobit.h"
#include <assert.h>
#include <limits.h>
#include <string.h>

struct assembler_t_
{
    /* kmer bit-mask used in assembly */
    kmer_t assemble_kmer_mask;

    /* kmer bit-mask used in alignment */
    kmer_t align_kmer_mask;

    /* read set */
    seqset_t* S;

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
};


assembler_t* assembler_alloc(size_t assemble_k, size_t align_k)
{
    assembler_t* A = malloc_or_die(sizeof(assembler_t));
    assert(A != NULL);

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


    // TODO: set n and m in some principled way
    A->B = bloom_alloc(8388608, 8);

    A->x = twobit_alloc();


    A->H = kmerhash_alloc();

    A->count_cutoff = 2;

    return A;
}


void assembler_free(assembler_t* A)
{
    seqset_free(A->S);
    bloom_free(A->B);
    kmerhash_free(A->H);
    twobit_free(A->x);
}



void assembler_add_seq(assembler_t* A, const char* seq, size_t seqlen)
{
    /* does the read contain non-nucleotide characters ? */
    size_t i;
    for (i = 0; i < seqlen; ++i) {
        /* XXX: for now, we just throw out any reads with non-nucleotide characters
         * */
        if (chartokmer(seq[i]) > 3) return;
    }

    // TODO: handle the case in which seqlen < k

    twobit_copy_n(A->x, seq, seqlen);
    seqset_inc(A->S, A->x);
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


static int seqset_value_cmp(const void* a_, const void* b_)
{
    seqset_value_t* a = (seqset_value_t*) a_;
    seqset_value_t* b = (seqset_value_t*) b_;

    if      (a->cnt < b->cnt) return 1;
    else if (a->cnt > b->cnt) return -1;
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
    /* An alignment is attempted from this many seeds spaced evenly across the
     * read. */
    const size_t seed_cnt = 3;

    /* We only consider at most the  first positions that a k-mer hashes to. */
    const size_t max_kmer_pos = 5;


    size_t i, j, seqlen;
    kmer_t x, y;

    kmer_pos_t* pos;
    size_t poslen;

    /* create an alignment structure for each contig */
    sw_t** sws = malloc_or_die(contigs_len * sizeof(sw_t*));
    for (i = 0; i < contigs_len; ++i) {
        sws[i] = sw_alloc(contigs[i]);
    }

    /* create an alignment structure for the reverse complement of each contig */
    twobit_t* rc_seq = twobit_alloc();
    sw_t** sws_rc = malloc_or_die(contigs_len * sizeof(sw_t*));
    for (i = 0; i < contigs_len; ++i) {
        twobit_revcomp(rc_seq, contigs[i]);
        sws_rc[i] = sw_alloc(rc_seq);
    }
    twobit_free(rc_seq);


    /* align every read! */

    sw_t* sw;
    int qpos, spos;


    /* minimum cost alignment */
    int aln_score, min_aln_score;
    sw_alignment_t aln;
    int aln_strand;

    memset(&aln, 0, sizeof(sw_alignment_t));

    size_t aln_cnt = 0;

    for (i = 0; i < xs_len; ++i) {
        seqlen = twobit_len(xs[i].seq);
        assert(seqlen >= A->align_k);

        aln.len = 0;
        min_aln_score = INT_MAX;

        for (j = 0; j < seed_cnt; ++j) {
            qpos = (int) (j * (seqlen - A->align_k) / (seed_cnt - 1));
            x = twobit_get_kmer(xs[i].seq, qpos, A->align_k);
            y = kmer_canonical(x, A->align_k);
            poslen = kmerhash_get(A->H, y, &pos);
            if (poslen > max_kmer_pos) poslen = max_kmer_pos;

            while (poslen--) {
                if (pos->contig_pos >= 0) {
                    if (x == y) {
                        spos = pos->contig_pos;
                        sw   = sws[pos->contig_idx];
                    }
                    else {
                        spos = (int) twobit_len(contigs[pos->contig_idx]) - pos->contig_pos - A->align_k;
                        sw   = sws_rc[pos->contig_idx];
                    }
                }
                else {
                    if (x == y) {
                        spos = (int) twobit_len(contigs[pos->contig_idx]) + pos->contig_pos - (A->align_k - 1);
                        sw   = sws_rc[pos->contig_idx];
                    }
                    else {
                        spos = -pos->contig_pos - 1;
                        sw   = sws[pos->contig_idx];
                    }
                }

                aln_score = sw_seeded_align(sw, xs[i].seq, spos, qpos, A->align_k);

                if (aln_score >= 0 && min_aln_score > aln_score) {
                    min_aln_score = aln_score;
                    sw_trace(sw, &aln);

                    if (pos->contig_pos >= 0) aln_strand = x == y ? 0 : 1;
                    else                      aln_strand = x == y ? 1 : 0;
                }

                pos++;
            }
        }

        // TODO: 
        // If there is a good alignment, write it.
        // Otherwise, compress the sequence.

        // XXX: report some rough statistics
        if (min_aln_score < (int) (seqlen / 2)) {
            ++aln_cnt;
        }
    }

    fprintf(stderr, "%zu / %zu (%0.2f%%) reads aligned.\n",
                    aln_cnt, xs_len, 100.0 * (double) aln_cnt / (double) xs_len);


    for (i = 0; i < contigs_len; ++i) {
        sw_free(sws[i]);
        sw_free(sws_rc[i]);
    }

    free(sws);
    free(sws_rc);

    free(aln.ops);
}



void assembler_assemble(assembler_t* A,
                        quip_block_writer_t writer,
                        void* writer_data)
{
    /* dump reads and sort by abundance */
    seqset_value_t* xs = seqset_dump(A->S);
    qsort(xs, seqset_size(A->S), sizeof(seqset_value_t), seqset_value_cmp);

    /* count kmers */
    size_t i, j, n, seqlen;
    kmer_t x, y;

    n = seqset_size(A->S);
    for (i = 0; i < n; ++i) {
        seqlen = twobit_len(xs[i].seq);
        x = 0;
        for (j = 0; j < seqlen; ++j) {
            x = ((x << 2) | twobit_get(xs[i].seq, j)) & A->assemble_kmer_mask;

            if (j + 1 >= A->assemble_k) {
                y = kmer_canonical(x, A->assemble_k);
                bloom_add(A->B, y, xs[i].cnt);
            }
        }
    }


    /* assemble contigs */
    size_t contigs_size = 512;
    size_t contigs_len  = 0;
    twobit_t** contigs = malloc_or_die(contigs_size * sizeof(twobit_t*));

    twobit_t* contig = twobit_alloc();

    for (i = 0; i < n && xs[i].cnt >= A->count_cutoff; ++i) {
        assembler_make_contig(A, xs[i].seq, contig);
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

    /* TODO: write compressed contigs */

    /* align reads to contigs */
    index_contigs(A, contigs, contigs_len);
    align_to_contigs(A, contigs, contigs_len, xs, n);


    for (i = 0; i < contigs_len; ++i) twobit_free(contigs[i]);
    free(contigs);
}


void assembler_flush(assembler_t* A)
{
    // TODO
}


