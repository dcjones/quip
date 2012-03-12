
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
#include <math.h>
#include <pthread.h>


/* Contigs are padded on both sides by this many As. */
static const size_t contig_padding = 50;

/* Maximum score for an alignment to be reported. */
static const double max_align_score = 1.20;

/* Number of alignment threads to use. */
static const size_t num_aln_threads = 2;


/* alignments */
typedef struct align_t_
{
    uint32_t contig_idx;
    uint8_t  strand;
    double   aln_score;
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

    /* an aligner for each contig*/
    sw_t** contig_aligners;
    sw_t** contig_rc_aligners;

    /* read being aligned */
    const twobit_t* query_seq;

    /* alignment threads */
    pthread_t*     aln_threads;
    pthread_cond_t aln_thread_cond;
    pthread_cond_t  seeds_aligned_cond;
    pthread_mutex_t seeds_aligned_mut;
    size_t threads_running;

    /* alignment seeds */
    kmer_pos_t* seeds;
    size_t seed_cnt;
    kmer_t seed_kmer, seed_kmer_c;
    int seed_qpos;

    /* next seed to try aligning the read to */
    size_t seed_idx;
    pthread_mutex_t seed_idx_mut;

    /* an alignment */
    align_t aln;
    pthread_mutex_t aln_mut;

    /* allocated size of the contigs array */
    size_t contigs_size;

    /* number of contigs stored in the contigs array */
    size_t contigs_len;

    /* should we try attempt to assemble contigs with the current batch of reads
     * */
    bool assemble_batch;

    /* nothing has been written yet */
    bool initial_state;
};


static void build_kmer_hash(assembler_t* A);
static void* aligner_thread(void*);


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
    A->initial_state = true;

    /* If we are not assembling, we do not need any of the data structure
     * initialized below. */
    if (!quick) {
        A->S = seqset_alloc();

        A->ord_size = 1024;
        A->ord = malloc_or_die(A->ord_size * sizeof(const twobit_t*));

        A->N = 0;

        A->B = bloom_alloc(2097152, 8);

        A->x = twobit_alloc();

        A->H = kmerhash_alloc();

        A->count_cutoff = 2;

        A->contigs_size = 512;
        A->contigs_len  = 0;
        A->contigs = malloc_or_die(A->contigs_size * sizeof(twobit_t*));

        A->contig_aligners    = NULL;
        A->contig_rc_aligners = NULL;

        A->aln.a.ops = NULL;

        A->query_seq = NULL;
        A->seeds = NULL;
        A->seed_cnt = 0;
        A->seed_idx = 0;

        pthread_mutex_init(&A->seed_idx_mut, NULL);
        pthread_mutex_init(&A->aln_mut, NULL);
        pthread_cond_init(&A->aln_thread_cond, NULL);

        pthread_cond_init(&A->seeds_aligned_cond, NULL);
        pthread_mutex_init(&A->seeds_aligned_mut, NULL);
        pthread_mutex_lock(&A->seeds_aligned_mut);

        A->aln_threads = malloc_or_die(num_aln_threads * sizeof(pthread_t));
        for (i = 0; i < num_aln_threads; ++i) {
            pthread_create(&A->aln_threads[i], NULL, aligner_thread, (void*) A);
        }
    }
    else {
        A->seqenc = seqenc_alloc_encoder(writer, writer_data);
    }

    return A;
}


void assembler_free(assembler_t* A)
{
    size_t i;
    if (!A->quick) {
        A->query_seq = NULL;
        pthread_cond_broadcast(&A->aln_thread_cond);

        void* val;
        for (i = 0; i < num_aln_threads; ++i) {
            pthread_join(A->aln_threads[i], &val);
        }
        free(A->aln_threads);

        pthread_mutex_destroy(&A->seed_idx_mut);
        pthread_mutex_destroy(&A->aln_mut);
        pthread_cond_destroy(&A->aln_thread_cond);

        pthread_cond_destroy(&A->seeds_aligned_cond);
        pthread_mutex_destroy(&A->seeds_aligned_mut);

    }

    seqset_free(A->S);
    free(A->ord);
    bloom_free(A->B);
    kmerhash_free(A->H);
    twobit_free(A->x);
    seqenc_free(A->seqenc);

    for (i = 0; i < A->contigs_len; ++i) {
        twobit_free(A->contigs[i]);
        sw_free(A->contig_aligners[i]);
        sw_free(A->contig_rc_aligners[i]);
    }
    free(A->contigs);
    free(A->contig_aligners);
    free(A->contig_rc_aligners);
    free(A->aln.a.ops);

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


static void* aligner_thread(void* arg)
{
    assembler_t* A = (assembler_t*) arg;

    int qlen, spos, slen;
    uint8_t strand;
    kmer_pos_t* seed;
    sw_t* aligner;
    double aln_score;

    pthread_mutex_t wait_mut;
    pthread_mutex_init(&wait_mut, NULL);
    pthread_mutex_lock(&wait_mut);

    size_t i;
    while (true) {
        pthread_mutex_lock(&A->seeds_aligned_mut);
        i = A->seed_idx;
        A->seed_idx++;

        if (i >= A->seed_cnt) {
            if (--A->threads_running == 0) pthread_cond_signal(&A->seeds_aligned_cond);
            pthread_cond_wait(&A->aln_thread_cond, &A->seeds_aligned_mut);
            pthread_mutex_unlock(&A->seeds_aligned_mut);

            if (A->query_seq == NULL) break;
            continue;
        }
        else {
            pthread_mutex_unlock(&A->seeds_aligned_mut);

            seed = &A->seeds[i];
            slen = twobit_len(A->contigs[seed->contig_idx]);
            qlen = twobit_len(A->query_seq);

            if (seed->contig_pos >= 0) {
                if (A->seed_kmer == A->seed_kmer_c) {
                    spos = seed->contig_pos;
                    strand = 0;
                }
                else {
                    spos = slen - seed->contig_pos - A->align_k;
                    strand = 1;
                }
            }
            else {
                if (A->seed_kmer == A->seed_kmer_c) {
                    spos = slen + seed->contig_pos - A->align_k + 1;
                    strand = 1;
                }
                else {
                    spos = -seed->contig_pos - 1;
                    strand = 0;
                }
            }

            /* Is a full alignment possible with this seed? */
            if ((strand == 0 && (spos < A->seed_qpos || slen - spos < qlen - A->seed_qpos)) ||
                (strand == 1 && (spos < qlen - A->seed_qpos || slen - spos < A->seed_qpos))) {
                continue;
            }

            if (strand == 0) {
                aligner = A->contig_aligners[seed->contig_idx];
            }
            else {
                aligner = A->contig_rc_aligners[seed->contig_idx];
            }

            sw_lock(aligner);

            aln_score = sw_seeded_align(
                    aligner,
                    A->query_seq,
                    spos,
                    A->seed_qpos,
                    A->align_k);

            pthread_mutex_lock(&A->aln_mut);
            if (aln_score >= 0 &&
                aln_score <= max_align_score &&
                aln_score < A->aln.aln_score)
            {
                A->aln.contig_idx = seed->contig_idx;
                A->aln.strand     = strand;
                A->aln.aln_score  = aln_score;
                sw_trace(aligner, &A->aln.a);
            }

            sw_unlock(aligner);
            pthread_mutex_unlock(&A->aln_mut);
        }
    }

    pthread_mutex_destroy(&wait_mut);
    return NULL;
}


static int align_read(assembler_t* A, const twobit_t* seq)
{
    int qlen = twobit_len(seq);
    A->query_seq = seq;
    A->aln.aln_score = HUGE_VAL;

    size_t i;
    for (i = 0; i < 3; ++i) {
        if      (i == 0) A->seed_qpos = 0;
        else if (i == 1) A->seed_qpos = (qlen - A->align_k) / 2;
        else if (i == 2) A->seed_qpos = qlen - A->align_k - 1;

        A->seed_kmer = twobit_get_kmer(
                            seq,
                            A->seed_qpos,
                            A->align_k);
        A->seed_kmer_c = kmer_canonical(A->seed_kmer, A->align_k);

        A->seed_cnt = kmerhash_get(A->H, A->seed_kmer_c, &A->seeds);
        A->seed_idx = 0;

        if (A->seed_cnt > 0) {

            // use this, if runtime is still an issue
            // poslen = poslen > max_seeds ? max_seeds : poslen;
     
            A->threads_running = num_aln_threads; 
            pthread_cond_broadcast(&A->aln_thread_cond);

            while (A->threads_running > 0) {
                pthread_cond_wait(&A->seeds_aligned_cond, &A->seeds_aligned_mut);
            }
        }
   }

    if (A->aln.aln_score < HUGE_VAL) {
        seqenc_encode_alignment(A->seqenc,
                A->aln.contig_idx, A->aln.strand, &A->aln.a, seq);

        return 1;
    }
    else {
        seqenc_encode_twobit_seq(A->seqenc, seq);
        return 0;
    }
}


/* Try desperately to align the given read. If no good alignment is found, just
 * encode the sequence. Return 0 if no alignment was found, 1 if a seed was
 * found but no good alignment, and 2 if a good alignment was found. */

 #if 0
static int align_read(assembler_t* A, const twobit_t* seq)
{
    /* We only consider the first few seed hits found in the hash table. The
     * should be in an approximately random order. */
    static const size_t max_seeds = 5;

    int ret = 0;

    /* position of the seed with the subject and query sequecne, resp. */
    int spos, qpos;

    /* subject and query lengths, reps */
    int slen, qlen;
    kmer_t x, y;
    uint8_t strand;

    /* positions matching the seed k-mer */
    kmer_pos_t* pos;
    size_t poslen;

    /* number of skipped seed candidates */
    size_t skipped;

    /* */
    double aln_score;
    A->aln.aln_score = HUGE_VAL;

    /* Don't try to align any reads that are shorer than the seed length */
    qlen = twobit_len(seq);
    if ((size_t) qlen < A->align_k) {
        seqenc_encode_twobit_seq(A->seqenc, seq);
        return false;
    }

    sw_t* aligner;

    size_t i, j;
    for (i = 0; i < 3; ++i) {
        if      (i == 0) qpos = 0;
        else if (i == 1) qpos = (qlen - A->align_k) / 2;
        else if (i == 2) qpos = qlen - A->align_k - 1;

        x = twobit_get_kmer(
                seq,
                qpos,
                A->align_k);
        y = kmer_canonical(x, A->align_k);

        poslen = kmerhash_get(A->H, y, &pos);
        poslen = poslen > max_seeds ? max_seeds : poslen;

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


            if (strand == 0) {
                aligner = A->contig_aligners[pos->contig_idx];
            }
            else {
                aligner = A->contig_rc_aligners[pos->contig_idx];
            }

            aln_score = sw_seeded_align(
                    aligner,
                    seq,
                    spos,
                    qpos,
                    A->align_k);

            if (aln_score >= 0 &&
                aln_score < A->aln.aln_score &&
                aln_score <= max_align_score)
            {
                A->aln.contig_idx = pos->contig_idx;
                A->aln.strand     = strand;
                A->aln.aln_score  = aln_score;
                sw_trace(aligner, &A->aln.a);
                ret = 2;
            }
            else if (aln_score >= 0 && ret == 0) ret = 1;
        }
    }

    if (A->aln.aln_score < INT_MAX) {
        seqenc_encode_alignment(A->seqenc,
                A->aln.contig_idx, A->aln.strand, &A->aln.a, seq);

        return ret;
    }
    else {
        seqenc_encode_twobit_seq(A->seqenc, seq);

        return ret;
    }
}
#endif


void assembler_add_seq(assembler_t* A, const char* seq, size_t seqlen)
{
    if (A->quick) {
        /* output the number of contigs (i.e., 0) */
        if (A->initial_state) {
            write_uint32(A->writer, A->writer_data, 0);
            A->initial_state = false;
        }

        seqenc_encode_char_seq(A->seqenc, seq, seqlen);
        return;
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
        /* output the number of contigs (i.e., 0) */
        if (A->initial_state) {
            write_uint32(A->writer, A->writer_data, 0);
            A->initial_state = false;
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
    kmer_t x = twobit_get_kmer(seed, 0, A->assemble_k);
    size_t i;
    for (i = A->assemble_k; i < twobit_len(seed); ++i) {
        bloom_ldec(A->B, kmer_canonical((x << 2) | twobit_get(seed, i), A->assemble_k));
    }


    /* expand the contig as far left as possible */
    unsigned int cnt, cnt2, cnt_best, cnt2_best;

    kmer_t nt, nt2, nt_best = 0, xc, y, z;


    /* Greedily append nucleotides to the contig using 
     * approximate k-mer counts stored in the bloom filter. */
    x = twobit_get_kmer(seed, 0, A->assemble_k);
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

    x = twobit_get_kmer(seed, twobit_len(seed) - A->assemble_k, A->assemble_k);
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

}


static void index_contigs(assembler_t* A)
{
    if (quip_verbose) fprintf(stderr, "indexing contigs ... ");

    build_kmer_hash(A);

    /* build aligner objects */
    A->contig_aligners    = malloc_or_die(A->contigs_len * sizeof(sw_t*));
    A->contig_rc_aligners = malloc_or_die(A->contigs_len * sizeof(sw_t*));

    size_t i;
    twobit_t* contig_rc = twobit_alloc();
    for (i = 0; i < A->contigs_len; ++i) {
        A->contig_aligners[i] = sw_alloc(A->contigs[i]);

        twobit_revcomp(contig_rc, A->contigs[i]);
        A->contig_rc_aligners[i] = sw_alloc(contig_rc);
    }
    twobit_free(contig_rc);

    if (quip_verbose) fprintf(stderr, "done.\n");
}

/* build new indexes after updating the consensus sequences */
static void reindex_contigs(assembler_t* A)
{
    if (quip_verbose) fprintf(stderr, "indexing contigs ... ");

    build_kmer_hash(A);

    size_t i;
    twobit_t* contig_rc = twobit_alloc();
    for (i = 0; i < A->contigs_len; ++i) {
        twobit_revcomp(contig_rc, A->contigs[i]);
        sw_set_subject(A->contig_aligners[i], A->contigs[i]);
        sw_set_subject(A->contig_rc_aligners[i], contig_rc);;
    }
    twobit_free(contig_rc);

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

    twobit_t* contig = twobit_alloc();
    twobit_t* padded_contig = twobit_alloc();
    kmer_t x, y;
    size_t len;

    size_t i, j;
    for (i = 0; i < n && xs[i].cnt >= A->count_cutoff; ++i) {
        if (!xs[i].is_twobit) continue;

        make_contig(A, xs[i].seq.tb, contig);

        /* reject over terribly short contigs */
        len = twobit_len(contig);
        if (len < twobit_len(xs[i].seq.tb) + 4 * A->assemble_k) {
            
            /* reclaim k-mers from the failed contig */
            x = 0;
            for (j = 0; j < len; ++j) {
                x = ((x << 2) | twobit_get(contig, j)) & A->assemble_kmer_mask;

                if (j + 1 >= A->assemble_k) {
                    y = kmer_canonical(x, A->assemble_k);
                    /* ideally we would add the k-mer back with its original
                     * count, but that information is lost. */
                    bloom_add(A->B, y, xs[i].cnt);
                }
            }

            continue;
        }

        twobit_clear(padded_contig);
        for (j = 0; j < contig_padding; ++j) {
            twobit_append_kmer(padded_contig, 0, 1);
        }

        twobit_append_twobit(padded_contig, contig);

        for (j = 0; j < contig_padding; ++j) {
            twobit_append_kmer(padded_contig, 0, 1);
        }

        if (A->contigs_len == A->contigs_size) {
            A->contigs_size *= 2;
            A->contigs = realloc_or_die(A->contigs, A->contigs_size * sizeof(twobit_t*));
        }

        A->contigs[A->contigs_len++] = twobit_dup(padded_contig);
    }

    twobit_free(contig);
    twobit_free(padded_contig);

    if (quip_verbose) fprintf(stderr, "done. (%zu contigs)\n", A->contigs_len);
}



void assembler_assemble(assembler_t* A)
{
    if (!A->quick && A->assemble_batch) {

        /* dump reads and sort by copy number*/
        seqset_value_t* xs = seqset_dump(A->S);
        qsort(xs, seqset_size(A->S), sizeof(seqset_value_t), seqset_value_cnt_cmp);

        size_t n = seqset_size(A->S);

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
        seqenc_set_contigs(A->seqenc, A->contigs, A->contigs_len);

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

        align_to_contigs(A, xs, n);

        free(xs);
        seqset_free(A->S);
        A->S = NULL;

        free(A->ord);
        A->ord = NULL;
    }

    seqenc_flush(A->seqenc);

    if (!A->quick) {
        seqenc_get_contig_consensus(A->seqenc, A->contigs);
        reindex_contigs(A);
    }

    A->initial_state = true;
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



