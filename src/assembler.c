
#include "assembler.h"
#include "bloom.h"
#include "hash.h"
#include "kmer.h"
#include "misc.h"
#include "seqset.h"
#include "twobit.h"
#include <assert.h>
#include <math.h>

struct assembler_t_
{
    /* kmer bit-mask */
    kmer_t mask;

    /* read set */
    seqset_t* S;

    /* current read */
    twobit_t* x;

    size_t k; /* k-mer size used for assembly */

    bloom_t*     B; /* k-mer table */

    /* count required to be nominated as a seed candidate */
    unsigned int count_cutoff;
};


assembler_t* assembler_alloc(size_t k)
{
    assembler_t* A = malloc_or_die(sizeof(assembler_t));
    assert(A != NULL);

    A->mask = 0;
    size_t i;
    for (i = 0; i < k; ++i) {
#ifdef WORDS_BIGENDIAN
        A->mask = (A->mask >> 2) | 0x3;
#else
        A->mask = (A->mask << 2) | 0x3;
#endif
    }

    A->S = seqset_alloc();

    A->k = k;
    // TODO: set n and m in some principled way
    A->B = bloom_alloc(8388608, 8);

    A->x = twobit_alloc();

    A->count_cutoff = 2;

    return A;
}


void assembler_free(assembler_t* A)
{
    seqset_free(A->S);
    bloom_free(A->B);
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

    // XXX
    /*twobit_append_twobit(contig, seed);*/
    /*return;*/

    /* expand the contig as far left as possible */
    unsigned int cnt, cnt_best = 0;

    kmer_t nt, nt_best = 0, x, xc, y;

    /* TODO: delete alle kmers in seed */

    x = twobit_get_kmer(seed, 0, A->k);
    while (true) {
        bloom_del(A->B, kmer_canonical(x, A->k));

#if WORDS_BIGENDIAN
        x = (x << 2) & A->mask;
#else 
        x = (x >> 2) & A->mask;
#endif
        cnt_best = 0;
        for (nt = 0; nt < 4; ++nt) {
#if WORDS_BIGENDIAN
            y = nt >> (2 * (A->k - 1));
#else
            y = nt << (2 * (A->k - 1));
#endif
            xc = kmer_canonical(x | y, A->k);
            cnt = bloom_get(A->B, xc);

            if (cnt > cnt_best) {
                cnt_best = cnt;
                nt_best  = nt;
            }
        }

        if (cnt_best > 0) {
#if WORDS_BIGENDIAN
            y = nt_best >> (2 * (A->k - 1));
#else
            y = nt_best << (2 * (A->k - 1));
#endif
            x = x | y;
            twobit_append_kmer(contig, nt_best, 1);
        }
        else break;
    }

    twobit_reverse(contig);
    twobit_append_twobit(contig, seed);

    x = twobit_get_kmer(seed, twobit_len(seed) - A->k, A->k);
    while (true) {
        bloom_del(A->B, kmer_canonical(x, A->k));

#if WORDS_BIGENDIAN
        x = (x >> 2) & A->mask;
#else 
        x = (x << 2) & A->mask;
#endif
        cnt_best = 0;
        for (nt = 0; nt < 4; ++nt) {
            xc = kmer_canonical(x | nt, A->k);
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


static void assembler_assemble(assembler_t* A)
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
#ifdef WORDS_BIGENDIAN
            x = ((x >> 2) | twobit_get(xs[i].seq, j)) & A->mask;
#else
            x = ((x << 2) | twobit_get(xs[i].seq, j)) & A->mask;
#endif

            if (j + 1 >= A->k) {
                y = kmer_canonical(x, A->k);
                bloom_add(A->B, y, xs[i].cnt);
            }
        }
    }


    fprintf(stdout, "%zu unique reads\n", seqset_size(A->S));

    FILE* f = fopen("contig.fa", "w");
    twobit_t* contig = twobit_alloc();

    for (i = 0; i < n && xs[i].cnt >= A->count_cutoff; ++i) {
        assembler_make_contig(A, xs[i].seq, contig);
        if (twobit_len(contig) <= A->k) continue;
        /* TODO: when we discard a contig, it would be nice if we could return
         * its k-mers to the bloom filter */

        fprintf(f, ">contig_%05zu\n", i);
        twobit_print(contig, f);
        fprintf(f, "\n\n");
    }

    twobit_free(contig);
    fclose(f);
}


void assembler_write(assembler_t* A, FILE* fout)
{
    assembler_assemble(A);

    fprintf(fout, "IOU: compressed data\n");
}



