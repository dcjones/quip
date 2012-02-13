
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



/* alignments seeds */
typedef struct aln_seed_t_
{
    uint32_t seq_idx;
    int      qpos;
    int      spos;
    uint8_t  strand;
    struct aln_seed_t_* next;
} aln_seed_t;


static void free_seeds(aln_seed_t** seeds, size_t n)
{
    aln_seed_t* seed;
    size_t i;
    for (i = 0; i < n; ++i) {
        seed = seeds[i];
        while (seed) {
            seed = seeds[i]->next;
            free(seeds[i]);
            seeds[i] = seed;
        }
    }
    free(seeds);
}


/* alignments */
typedef struct align_t_
{
    uint32_t contig_idx;
    uint8_t  strand;
    int      aln_score;
    sw_alignment_t a;
} align_t;



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

    /* assembled contigs */
    twobit_t** contigs;

    /* allocated size of the contigs array */
    size_t contigs_size;

    /* number of contigs stored in the contigs arary */
    size_t contigs_len;

    /* should we try attempt to assemble contigs with the current batch of reads
     * */
    bool assemble_batch;

    /* nothing has been written yet */
    bool initial_state;
};


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

        // TODO: set n and m in some principled way
        A->B = bloom_alloc(4194304, 8);

        A->x = twobit_alloc();

        A->H = kmerhash_alloc();

        A->count_cutoff = 2;

        A->contigs_size = 512;
        A->contigs_len = 0;;
        A->contigs = malloc_or_die(A->contigs_size * sizeof(twobit_t*));
    }
    else {
        A->seqenc = seqenc_alloc_encoder(writer, writer_data);
        A->initial_state = true;
    }

    return A;
}


static void assembler_reset(assembler_t* A)
{
    seqset_clear(A->S);
    memset(A->ord, 0, A->ord_size * sizeof(uint32_t));
    A->N = 0;

    if (A->B) bloom_clear(A->B);
}


void assembler_free(assembler_t* A)
{
    seqset_free(A->S);
    free(A->ord);
    bloom_free(A->B);
    kmerhash_free(A->H);
    twobit_free(A->x);
    seqenc_free(A->seqenc);

    size_t i;
    for (i = 0; i < A->contigs_len; ++i) {
        twobit_free(A->contigs[i]);
    }
    free(A->contigs);

    free(A);
}


void assembler_clear_contigs(assembler_t* A)
{
    size_t i;
    for (i = 0; i < A->contigs_len; ++i) {
        twobit_free(A->contigs[i]);
    }
    A->contigs_len = 0;
}



void assembler_add_seq(assembler_t* A, const char* seq, size_t seqlen)
{
    if (A->quick) {
        if (A->initial_state) {
            /* output the number of contigs (i.e., 0) */
            write_uint32(A->writer, A->writer_data, 0);
            A->initial_state = false;
        }

        seqenc_encode_char_seq(A->seqenc, seq, seqlen);
        return;
    }

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



static void make_contig(assembler_t* A, twobit_t* seed, twobit_t* contig)
{
    twobit_clear(contig);


    /* delete all kmers in the seed */
    kmer_t x = twobit_get_kmer(seed, 0, A->assemble_k);
    size_t i;
    for (i = A->assemble_k; i < twobit_len(seed); ++i) {
        bloom_dec(A->B, kmer_canonical((x << 2) | twobit_get(seed, i), A->assemble_k));
    }


    /* expand the contig as far left as possible */
    unsigned int cnt, cnt_best = 0;

    kmer_t nt, nt_best = 0, xc, y;


    x = twobit_get_kmer(seed, 0, A->assemble_k);
    while (true) {
        bloom_dec(A->B, kmer_canonical(x, A->assemble_k));

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
        bloom_dec(A->B, kmer_canonical(x, A->assemble_k));

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


static void index_contigs(assembler_t* A)
{
    fprintf(stderr, "indexing contigs ... ");
    size_t i;
    size_t len;
    size_t pos;
    kmer_t x, y;
    twobit_t* contig;
    for (i = 0; i < A->contigs_len; ++i) {
        contig = A->contigs[i];
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



static aln_seed_t** seed_alignments(assembler_t* A,
                            seqset_value_t* xs, size_t xs_len)
{
    if (verbose) fprintf(stderr, "\tseeding alignments ... ");

    aln_seed_t* seed;
    aln_seed_t** seeds= malloc_or_die(A->contigs_len * sizeof(aln_seed_t*));
    memset(seeds, 0, A->contigs_len * sizeof(aln_seed_t*));


    /* We only consider the first few seed hits found in the hash table. The
     * should be in an approximately random order. */
    const size_t max_seeds = 5;

    size_t poslen, skipped;
    kmer_pos_t* pos;
    kmer_t x, y;
    size_t i, j, k;

    int spos, qpos;
    int slen, qlen;
    uint8_t strand;


    for (i = 0; i < xs_len; ++i) {

        /* Currently, we do not bother aligning reads with Ns (that cannot
         * be encoded in twobit). */
        if (!xs[i].is_twobit) continue;

        /* Don't try to align any reads that are shorer than the seed length */
        qlen = twobit_len(xs[i].seq.tb);
        if ((size_t) qlen < A->align_k) continue;

        /* loop over seed positions */
        for (k = 0; k < 3; ++k) {
            if      (k == 0) qpos = 0;
            else if (k == 1) qpos = (qlen - A->align_k) / 2;
            else if (k == 2) qpos = qlen - A->align_k - 1;

            x = twobit_get_kmer(
                    xs[i].seq.tb,
                    qpos,
                    A->align_k);
            y = kmer_canonical(x, A->align_k);

            poslen = kmerhash_get(A->H, y, &pos);
            poslen = poslen > max_seeds ? max_seeds : poslen;


            /* Adjust positions according to strand and record hits.  */
            for (j = 0, skipped = 0; j < poslen && j - skipped < max_seeds; ++j, ++pos) {
                slen = twobit_len(A->contigs[pos->contig_idx]);

                if (pos->contig_pos >= 0) {
                    if (x == y) {
                        spos = pos->contig_pos;
                        strand = 0;
                    }
                    else {
                        spos = slen - pos->contig_pos - A->align_k;
                        strand = 1;
                    }
                }
                else {
                    if (x == y) {
                        spos = slen + pos->contig_pos - A->align_k + 1;
                        strand = 1;
                    }
                    else {
                        spos = -pos->contig_pos - 1;
                        strand = 0;
                    }
                }

                /* Is a full alignment possible with this seed? */
                if ((strand == 0 && (spos < qpos || slen - spos < qlen - qpos)) ||
                    (strand == 1 && (spos < qlen - qpos || slen - spos < qpos))) {
                    ++skipped;
                    continue;
                }

                seed = malloc_or_die(sizeof(aln_seed_t));
                seed->spos    = spos;
                seed->qpos    = qpos;
                seed->strand  = strand;
                seed->seq_idx = i;
                seed->next    = seeds[pos->contig_idx];
                seeds[pos->contig_idx] = seed;
            }
        }
    }


    if (verbose) fprintf(stderr, "done.\n");

    return seeds;
}



static align_t* extend_seeds(assembler_t* A,
                             seqset_value_t* xs, size_t xs_len,
                             aln_seed_t** seeds)
{
    if (verbose) fprintf(stderr, "\tlocal alignment ... ");

    size_t i, j;
    aln_seed_t* seed;

    align_t* alns = malloc_or_die(xs_len * sizeof(align_t));
    for (i = 0; i < xs_len; ++i) {
        alns[i].aln_score = INT_MAX;
        memset(&alns[i].a, 0, sizeof(sw_alignment_t));
    }

    sw_t* sw;
    sw_t* sw_rc;
    twobit_t* contig_rc = twobit_alloc();
    size_t seqlen;

    int aln_score;

    for (i = 0, j = 0; i < A->contigs_len; ++i) {

        if (!seeds[i]) continue;

        sw = sw_alloc(A->contigs[i]);

        twobit_revcomp(contig_rc, A->contigs[i]);
        sw_rc = sw_alloc(contig_rc);

        while (seeds[i]) {

            seed = seeds[i];
            seqlen = twobit_len(xs[seed->seq_idx].seq.tb);

            /* if we already have a good alignment for this read, skip over it
             * */
            if (alns[seed->seq_idx].aln_score < 12 * (int) seqlen / 10) {
                seeds[i] = seed->next;
                free(seed);
                continue;
            }


            /* local alignment */
            aln_score = sw_seeded_align(
                            seed->strand == 0 ? sw : sw_rc,
                            xs[seed->seq_idx].seq.tb,
                            seed->spos,
                            seed->qpos,
                            A->align_k);


            /* better than the current alignment ? */
            if (aln_score >= 0 &&
                aln_score < 13 * (int) seqlen / 10 && /* crude cutoff for what an acceptable alignment is */
                aln_score < alns[seed->seq_idx].aln_score)
            {
                alns[seed->seq_idx].contig_idx = i;
                alns[seed->seq_idx].strand = seed->strand;
                alns[seed->seq_idx].aln_score = aln_score;
                sw_trace(seed->strand == 0 ? sw : sw_rc, &alns[seed->seq_idx].a);
            }

            seeds[i] = seed->next;
            free(seed);
        }

        sw_free(sw);
        sw_free(sw_rc);
    }

    twobit_free(contig_rc);

    if (verbose) fprintf(stderr, "done.");

    return alns;
}



static void align_to_contigs(assembler_t* A,
                             seqset_value_t* xs, size_t xs_len)
{
    if (verbose) fprintf(stderr, "aligning reads to contigs ...\n");

    /* Sort reads by index. Afterwards xs[i].idx == i should be true for every
     * i, where 0 <= i < xs_len */
    qsort(xs, xs_len, sizeof(seqset_value_t), seqset_value_idx_cmp);

    /* First pass: hash k-mers against the contigs to generate alignment seeds */
    aln_seed_t** seeds = seed_alignments(A, xs, xs_len);

    /* Second pass: one contig at a time, perform local alignment. */
    align_t* alns = extend_seeds(A, xs, xs_len, seeds);

    /* Lastly: output compressed reads */
    size_t aln_count = 0;
    const char* ebseq;
    size_t i;
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
                ebseq = xs[A->ord[i]].seq.eb;
                seqenc_encode_char_seq(A->seqenc, ebseq, strlen(ebseq));
            }
        }
    }

    if (verbose) {
        fprintf(stderr,
                "\t%zu / %zu [%0.2f%%] reads aligned to contigs\n",
                aln_count, (size_t) A->N, 100.0 * (double) aln_count / (double) A->N);
    }


    free_seeds(seeds, A->contigs_len);

    for (i = 0; i < xs_len; ++i) {
        free(alns[i].a.ops);
    }
    free(alns);

    if (verbose) fprintf(stderr, "done.\n");
}



/* Count the number of occurences of each k-mer in the array of reads xs, of
 * length n. */
static void count_kmers(assembler_t* A, seqset_value_t* xs, size_t n)
{
    if (verbose) fprintf(stderr, "counting k-mers ... ");

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

    if (verbose) fprintf(stderr, "done.\n");
}


/* Heuristically build contigs from k-mers counts and a set of reads */
static void make_contigs(assembler_t* A, seqset_value_t* xs, size_t n)
{
    if (verbose) fprintf(stderr, "assembling contigs ... ");

    twobit_t* contig = twobit_alloc();
    kmer_t x, y;
    size_t len;

    size_t i, j;
    for (i = 0; i < n && xs[i].cnt >= A->count_cutoff; ++i) {
        if (!xs[i].is_twobit) continue;

        make_contig(A, xs[i].seq.tb, contig);

        /* skip overy terribly short contigs */
        len = twobit_len(contig);
        if (len < 3 * A->assemble_k) {
            
            /* reclaim k-mers from the failed contig */
            x = 0;
            for (j = 0; j < len; ++j) {
                x = ((x << 2) | twobit_get(contig, j)) & A->assemble_kmer_mask;

                if (j + 1 >= A->assemble_k) {
                    y = kmer_canonical(x, A->assemble_k);
                    /* ideally we would add the k-mer back with its original
                     * count, but that information is lost. */
                    bloom_add(A->B, y, 1);
                }
            }

            continue;
        }


        if (A->contigs_len == A->contigs_size) {
            A->contigs_size *= 2;
            A->contigs = realloc_or_die(A->contigs, A->contigs_size * sizeof(twobit_t*));
        }

        A->contigs[A->contigs_len++] = twobit_dup(contig);
    }

    twobit_free(contig);

    if (verbose) fprintf(stderr, "done. (%zu contigs)\n", A->contigs_len);
}



void assembler_assemble(assembler_t* A)
{
    if (A->quick) {
        seqenc_flush(A->seqenc);
        A->initial_state = true;
        return;
    }

    /* dump reads and sort by copy number*/
    seqset_value_t* xs = seqset_dump(A->S);
    qsort(xs, seqset_size(A->S), sizeof(seqset_value_t), seqset_value_cnt_cmp);

    size_t n = seqset_size(A->S);

    if (A->assemble_batch) {

        count_kmers(A, xs, n);
        make_contigs(A, xs, n);

        /* Only assemble the first batch. */
        A->assemble_batch = false;

        /* And free up memory we won't need any more */
        bloom_free(A->B);
        A->B = NULL;

        /* Now that we have some memory to spare, bring the sequence encoder
         * online. */
        A->seqenc = seqenc_alloc_encoder(A->writer, A->writer_data);

        /* write the number of contigs and their lengths  */
        size_t i;
        write_uint32(A->writer, A->writer_data, A->contigs_len);
        for (i = 0; i < A->contigs_len; ++i) {
            write_uint32(A->writer, A->writer_data, twobit_len(A->contigs[i]));
        }

        /* write the contigs */
        for (i = 0; i < A->contigs_len; ++i) {
            seqenc_encode_twobit_seq(A->seqenc, A->contigs[i]);
        }
            
        index_contigs(A);

    }
    else {
        /* indicate that no new contigs will be used */
        write_uint32(A->writer, A->writer_data, 0); 
    }

    if (!A->seqenc) {
        A->seqenc = seqenc_alloc_encoder(A->writer, A->writer_data);
    }


    align_to_contigs(A, xs, n);

    free(xs);

    seqenc_flush(A->seqenc);

    assembler_reset(A);
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
        uint32_t contig_count = read_uint32(D->reader, D->reader_data);
        uint32_t* contig_lens = malloc_or_die(contig_count * sizeof(uint32_t));
        uint32_t i;
        for (i = 0; i < contig_count; ++i) {
            contig_lens[i] = read_uint32(D->reader, D->reader_data);
        }

        if (contig_count > 0) {
            seqenc_prepare_decoder(D->seqenc, contig_count, contig_lens);
        }

        free(contig_lens);

        D->init_state = false;
    }

    seqenc_decode(D->seqenc, x, n);
}


void disassembler_reset(disassembler_t* D)
{
    seqenc_reset_decoder(D->seqenc);
    D->init_state = true;
}



