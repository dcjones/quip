
#include "assembler.h"
#include "bloom.h"
#include "hash.h"
#include "kmer.h"
#include "misc.h"
#include "twobit.h"
#include <assert.h>
#include <math.h>

struct assembler_t_
{
    /* kmer bit-mask */
    kmer_t mask;

    /* nucleotide sequences */
    twobit_t* xs;

    /* the sequence of read i is at xs[offset[i]] */
    uint32_t* offsets;
    uint32_t next_offset;

    /* size fo offsets and read_lens arrays */
    size_t max_n;

    size_t n; /* number of stored reads */
    size_t k; /* k-mer size used for assembly */

    bloom_t*     B; /* k-mer table */
    kmer_hash_t* H; /* seed candidates */

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

    A->k = k;
    // TODO: set n and m in some principled way
    A->B = bloom_alloc(4194304, 8);

    A->xs = twobit_alloc();

    A->n           = 0;
    A->max_n       = 512;
    A->offsets     = malloc_or_die(A->max_n * sizeof(uint32_t));
    A->offsets[0] = 0;

    A->H = kmer_hash_alloc();
    A->count_cutoff = 10;

    return A;
}


void assembler_free(assembler_t* A)
{
    bloom_free(A->B);
    twobit_free(A->xs);
    free(A->offsets);
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


    /* append to xs */
    while (A->n + 1 >= A->max_n) {
        A->max_n *= 2;
        A->offsets   = realloc_or_die(A->offsets, A->max_n * sizeof(uint32_t));
    }

    twobit_append_n(A->xs, seq, seqlen);
    A->offsets[A->n + 1] = A->offsets[A->n] + seqlen;
    ++A->n;



    /* hash kmers */
    size_t j;
    unsigned int cnt;
    kmer_t y, x = 0;
    for (i = A->offsets[A->n - 1], j = 0; i < A->offsets[A->n]; ++i, ++j) {
#ifdef WORDS_BIGENDIAN
        x = ((x >> 2) | twobit_get(A->xs, i)) & A->mask;
#else
        x = ((x << 2) | twobit_get(A->xs, i)) & A->mask;
#endif

        if (j + 1 >= A->k) {
            y = kmer_canonical(x, A->k);
            cnt = bloom_inc(A->B, y);
            if (cnt >= A->count_cutoff) {
                kmer_hash_put(A->H, y, cnt);
            }
        }
    }
}


typedef struct seed_t_
{
    kmer_t x;
    unsigned int count;
} seed_t;



static void assembler_make_contig(assembler_t* A, kmer_t seed, twobit_t* contig)
{
    twobit_clear(contig);

    /* expand the contig as far left as possible */
    unsigned int cnt0, cnt, cnt_best = 0;
    int min_diff, diff;

    cnt0 = bloom_get(A->B, kmer_canonical(seed, A->k));
    kmer_t nt, nt_best = 0, x, y;

    x = seed;
    while (true) {
        bloom_del(A->B, kmer_canonical(x, A->k));

#if WORDS_BIGENDIAN
        x = (x << 2) & A->mask;
#else 
        x = (x >> 2) & A->mask;
#endif
        min_diff = 1000000;
        for (nt = 0; nt < 4; ++nt) {
#if WORDS_BIGENDIAN
            y = nt >> (2 * (A->k - 1));
#else
            y = nt << (2 * (A->k - 1));
#endif
            cnt = bloom_get(A->B, kmer_canonical(x | y, A->k));
            diff = abs((int) cnt - (int) cnt0);
            if (diff < min_diff) {
                min_diff = diff;
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
    twobit_append_kmer(contig, seed, A->k);

    x = seed;
    while (true) {
        bloom_del(A->B, kmer_canonical(x, A->k));

#if WORDS_BIGENDIAN
        x = (x >> 2) & A->mask;
#else 
        x = (x << 2) & A->mask;
#endif
        min_diff = 1000000;
        cnt_best = 0;
        for (nt = 0; nt < 4; ++nt) {
            cnt = bloom_get(A->B, kmer_canonical(x | nt, A->k));
            diff = abs((int) cnt - (int) cnt0);
            if (diff < min_diff) {
                min_diff = diff;
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


static void assembler_assemble(assembler_t* A)
{
    kmer_count_pair_t* seeds = kmer_hash_dump_sorted(A->H);
    size_t num_seeds = kmer_hash_size(A->H);

    fprintf(stdout, "%zu seeds\n", num_seeds);

    FILE* f = fopen("contig.fa", "w");

    twobit_t* contig = twobit_alloc();
    size_t i;
    for (i = 0; i < num_seeds; ++i) {
        assembler_make_contig(A, seeds[i].x, contig);
        if (twobit_len(contig) <= A->k) continue;
        fprintf(f, ">contig_%05zu\n", i);
        twobit_print(contig, f);
        fprintf(f, "\n\n");
    }
    twobit_free(contig);

    fclose(f);

    free(seeds);
}


void assembler_write(assembler_t* A, FILE* fout)
{
    assembler_assemble(A);

    fprintf(fout, "IOU: compressed data\n");
}



