
#include "quipfmt.h"
#include "assembler.h"
#include "misc.h"
#include "qualenc.h"
#include "idenc.h"
#include "crc64.h"
#include <stdint.h>
#include <string.h>
#include <pthread.h>
#include <errno.h>

static const uint8_t quip_header_magic[6] = 
    {0xff, 'Q', 'U', 'I', 'P', 0x00};

static const uint8_t quip_header_version = 0x01;

static const size_t assembler_k = 25;
static const size_t aligner_k   = 12;

/* maximum number of bases per block */
static const size_t block_size = 65000000;

/* highest quality score supported. */
static const char qual_scale_size = 64;

/* Maximum number of sequences to read before they are compressed. */
#define chunk_size 5000 


static void write_uint8(quip_writer_t writer, void* writer_data, uint8_t x)
{
    writer(writer_data, &x, 1);
}


static void write_uint32(quip_writer_t writer, void* writer_data, uint32_t x)
{
    uint8_t bytes[4] = { (uint8_t) (x >> 24),
                         (uint8_t) (x >> 16),
                         (uint8_t) (x >> 8),
                         (uint8_t) x };
    writer(writer_data, bytes, 4);
}


static void write_uint64(quip_writer_t writer, void* writer_data, uint64_t x)
{
    uint8_t bytes[8] = { (uint8_t) (x >> 56),
                         (uint8_t) (x >> 48),
                         (uint8_t) (x >> 40),
                         (uint8_t) (x >> 32),
                         (uint8_t) (x >> 24),
                         (uint8_t) (x >> 16),
                         (uint8_t) (x >> 8),
                         (uint8_t) x };

    writer(writer_data, bytes, 8);
}


static uint8_t read_uint8(quip_reader_t reader, void* reader_data)
{
    uint8_t x;
    size_t cnt = reader(reader_data, &x, 1);

    if (cnt == 0) {
        fprintf(stderr, "Unexpected end of file.\n");
        exit(EXIT_FAILURE);
    }

    return x;
}


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


static uint64_t read_uint64(quip_reader_t reader, void* reader_data)
{
    uint8_t bytes[8];
    size_t cnt = reader(reader_data, bytes, 8);

    if (cnt < 8) {
        fprintf(stderr, "Unexpected end of file.\n");
        exit(EXIT_FAILURE);
    }

    return ((uint64_t) bytes[0] << 56) |
           ((uint64_t) bytes[1] << 48) |
           ((uint64_t) bytes[2] << 40) |
           ((uint64_t) bytes[3] << 32) |
           ((uint64_t) bytes[4] << 24) |
           ((uint64_t) bytes[5] << 16) |
           ((uint64_t) bytes[6] << 8) |
           ((uint64_t) bytes[7]);
}


static void pthread_join_or_die(pthread_t thread, void** value_ptr)
{
    int ret = pthread_join(thread, value_ptr);
    if (ret != 0) {
        fputs("pthread_join error: ", stderr);
        switch (ret)
        {
            case EDEADLK:
                fputs("deadlock detected.\n", stderr);
                break;

            case EINVAL:
                fputs("non-joinable thread joined.\n", stderr);
                break;

            case ESRCH:
                fputs("thread could not be found.\n", stderr);
                break;

            default:
                fputs("mysterious non-specific error.\n", stderr);
        }

        exit(EXIT_FAILURE);
    }
}

static void pthread_create_or_die(
    pthread_t* thread, const pthread_attr_t* attr, void *(*start_routine)(void*), void* arg)
{
    int ret = pthread_create(thread, attr, start_routine, arg);
    if (ret != 0) {
        fputs("pthread_create error: ", stderr);
        switch (ret) 
        {
            case EAGAIN:
                fputs("insufficient resources.\n", stderr);
                break;

            case EINVAL:
                fputs("invalid attributes.\n", stderr);
                break;

            default:
                fputs("mysterious non-specific error.\n", stderr);
        }

        exit(EXIT_FAILURE);
    }
}



struct quip_out_t_
{
    /* sequence buffers */
    short_read_t* chunk[chunk_size];
    size_t chunk_len;

    /* function for writing compressed data */
    quip_writer_t writer;
    void* writer_data;

    /* algorithms to compress ids, qualities, and sequences, resp. */
    idenc_t*     idenc;
    qualenc_t*   qualenc;
    assembler_t* assembler;

    /* number of bases encoded in the current block */
    size_t buffered_reads;
    size_t buffered_bases;

    /* block specific checksums */
    uint64_t id_crc;
    uint64_t seq_crc;
    uint64_t qual_crc;

    /* for compression statistics */
    uint32_t id_bytes;
    uint32_t qual_bytes;
    uint32_t seq_bytes;

    /* run length encoded read lengths */
    uint32_t* readlen_vals;
    uint32_t* readlen_lens;
    size_t readlen_count, readlen_size;

    /* run length encoded quality score scheme guesses */
    char*     qual_scheme_vals;
    uint32_t* qual_scheme_lens;
    size_t qual_scheme_count, qual_scheme_size;

    /* general statistics */
    uint64_t total_reads;
    uint64_t total_bases;

    /* Have all the reads been written? */
    bool finished;
};


static void* id_compressor_thread(void* ctx)
{
    quip_out_t* C = (quip_out_t*) ctx;

    size_t i;
    for (i = 0; i < C->chunk_len; ++i) {
        idenc_encode(C->idenc, C->chunk[i]);
        C->id_crc = crc64_update(
            (uint8_t*) C->chunk[i]->id1.s,
            C->chunk[i]->id1.n, C->id_crc);
    }

    return NULL;
}

static void* seq_compressor_thread(void* ctx)
{
    quip_out_t* C = (quip_out_t*) ctx;

    size_t i;
    for (i = 0; i < C->chunk_len; ++i) {
        assembler_add_seq(C->assembler, C->chunk[i]->seq.s, C->chunk[i]->seq.n);
        C->seq_crc = crc64_update(
            (uint8_t*) C->chunk[i]->seq.s,
            C->chunk[i]->seq.n, C->seq_crc);
    }

    return NULL;
}

static void* qual_compressor_thread(void* ctx)
{
    quip_out_t* C = (quip_out_t*) ctx;

    size_t i;
    for (i = 0; i < C->chunk_len; ++i) {
        qualenc_encode(C->qualenc, C->chunk[i]);
        C->qual_crc = crc64_update(
            (uint8_t*) C->chunk[i]->qual.s,
            C->chunk[i]->qual.n, C->qual_crc);
    }

    return NULL;
}


quip_out_t* quip_out_alloc(quip_writer_t writer, void* writer_data, bool quick)
{
    quip_out_t* C = malloc_or_die(sizeof(quip_out_t));
    C->writer = writer;
    C->writer_data = writer_data;

    size_t i;
    for (i = 0; i < chunk_size; ++i) {
        C->chunk[i] = fastq_alloc_seq();
    }
    C->chunk_len = 0;

    C->buffered_reads = 0;
    C->buffered_bases = 0;

    C->id_crc   = 0;
    C->seq_crc  = 0;
    C->qual_crc = 0;

    C->id_bytes   = 0;
    C->qual_bytes = 0;
    C->seq_bytes  = 0;

    C->readlen_size  = 1;
    C->readlen_count = 0;
    C->readlen_vals = malloc_or_die(C->readlen_size * sizeof(uint32_t));
    C->readlen_lens = malloc_or_die(C->readlen_size * sizeof(uint32_t));

    C->qual_scheme_size  = 1;
    C->qual_scheme_count = 1;
    C->qual_scheme_vals = malloc_or_die(C->qual_scheme_size * sizeof(char));
    C->qual_scheme_lens = malloc_or_die(C->qual_scheme_size * sizeof(uint32_t));
    C->qual_scheme_vals[0] = '!';

    C->total_reads = 0;
    C->total_bases = 0;

    C->idenc     = idenc_alloc_encoder(writer, (void*) writer_data);
    C->qualenc   = qualenc_alloc_encoder(writer, (void*) writer_data);
    C->assembler = assembler_alloc(writer, (void*) writer_data, assembler_k, aligner_k, quick);

    /* write header */
    C->writer(C->writer_data, quip_header_magic, 6);
    C->writer(C->writer_data, &quip_header_version, 1);

    C->finished = false;

    return C;
}

static void quip_out_add_readlen(quip_out_t* C, size_t l)
{
    if (C->readlen_count == 0 || C->readlen_vals[C->readlen_count - 1] != l) {
        if (C->readlen_count >= C->readlen_size) {
            while (C->readlen_count >= C->readlen_size) {
                C->readlen_size *= 2;
            }

            C->readlen_vals = realloc_or_die(C->readlen_vals, C->readlen_size * sizeof(uint32_t));
            C->readlen_lens = realloc_or_die(C->readlen_vals, C->readlen_size * sizeof(uint32_t));
        }

        C->readlen_vals[C->readlen_count] = l;
        C->readlen_lens[C->readlen_count] = 1;
        C->readlen_count++;
    }
    else C->readlen_lens[C->readlen_count - 1]++;
}


void quip_out_flush_block(quip_out_t* C)
{
    if (quip_verbose) {
        fprintf(stderr, "writing a block of %zu compressed bases...\n", C->buffered_bases);
    }

    /* write the number of encoded reads and bases in the block */
    write_uint32(C->writer, C->writer_data, C->buffered_reads);
    write_uint32(C->writer, C->writer_data, C->buffered_bases);

    /* write the run length encoded read lengths. */
    size_t i;
    for (i = 0; i < C->readlen_count; ++i) {
        write_uint32(C->writer, C->writer_data, C->readlen_vals[i]);
        write_uint32(C->writer, C->writer_data, C->readlen_lens[i]);
    }

    /* write the run length encoded quality scheme guesses */
    for (i = 0; i < C->qual_scheme_count; ++i) {
        write_uint8(C->writer, C->writer_data, C->qual_scheme_vals[i]);
        write_uint32(C->writer, C->writer_data, C->qual_scheme_lens[i]);
    }


    /* finish coding */
    size_t comp_id_bytes   = idenc_finish(C->idenc);
    size_t comp_seq_bytes  = assembler_finish(C->assembler);
    size_t comp_qual_bytes = qualenc_finish(C->qualenc);

    /* write the size of each chunk and its checksum */
    write_uint32(C->writer, C->writer_data, C->id_bytes);
    write_uint32(C->writer, C->writer_data, comp_id_bytes);
    write_uint64(C->writer, C->writer_data, C->id_crc);

    write_uint32(C->writer, C->writer_data, C->seq_bytes);
    write_uint32(C->writer, C->writer_data, comp_seq_bytes);
    write_uint64(C->writer, C->writer_data, C->seq_crc);

    write_uint32(C->writer, C->writer_data, C->qual_bytes);
    write_uint32(C->writer, C->writer_data, comp_qual_bytes);
    write_uint64(C->writer, C->writer_data, C->qual_crc);

    /* write compressed ids */
    idenc_flush(C->idenc);
    if (quip_verbose) {
        fprintf(stderr, "\tid: %u / %u (%0.2f%%)\n",
                (unsigned int) comp_id_bytes, (unsigned int) C->id_bytes,
                100.0 * (double) comp_id_bytes / (double) C->id_bytes);
    }

    /* write compressed sequences */
    assembler_flush(C->assembler);
    if (quip_verbose) {
        fprintf(stderr, "\tseq: %u / %u (%0.2f%%)\n",
                (unsigned int) comp_seq_bytes, (unsigned int) C->seq_bytes,
                100.0 * (double) comp_seq_bytes / (double) C->seq_bytes);
    }


    /* write compressed qualities */
    qualenc_flush(C->qualenc);
    if (quip_verbose) {
        fprintf(stderr, "\tqual: %u / %u (%0.2f%%)\n",
                (unsigned int) comp_qual_bytes, (unsigned int) C->qual_bytes,
                100.0 * (double) comp_qual_bytes / (double) C->qual_bytes);
    }

    /* reset */
    C->buffered_reads = 0;
    C->buffered_bases = 0;
    C->id_bytes       = 0;
    C->qual_bytes     = 0;
    C->seq_bytes      = 0;
    C->id_crc         = 0;
    C->seq_crc        = 0;
    C->qual_crc       = 0;
    C->readlen_count  = 0;

    C->qual_scheme_vals[0] = C->qual_scheme_vals[C->qual_scheme_count - 1];
    C->qual_scheme_lens[0] = 0;
    C->qual_scheme_count = 1;
}


/* Ensure the proper quality score scheme is used for the
 * current chunk. */
static void update_qual_scheme_guess(quip_out_t* C)
{
    size_t i, j;
    char last_base_qual = C->qual_scheme_vals[C->qual_scheme_count - 1];
    char min_qual = '~';
    char max_qual = '!';
    for (i = 0; i < C->chunk_len; ++i) {
        for (j = 0; j < C->chunk[i]->qual.n; ++j) {
            if (C->chunk[i]->qual.s[j] < min_qual) {
                min_qual = C->chunk[i]->qual.s[j];
            }

            if (C->chunk[i]->qual.s[j] > max_qual) {
                max_qual = C->chunk[i]->qual.s[j];
            }
        }
    }

    if (max_qual - min_qual > max_qual) {
        fprintf(stderr, "Invalid quality score scheme: are large range is used than quip "
                        "currently supports..\n");
        exit(EXIT_FAILURE);
    }

    if (min_qual <  last_base_qual ||
        max_qual >= last_base_qual + qual_scale_size) {

        /* new quality scheme guess */
        C->qual_scheme_count++;
        if (C->qual_scheme_count >= C->qual_scheme_size) {
            C->qual_scheme_size++;
            C->qual_scheme_vals =
                realloc_or_die(
                    C->qual_scheme_vals,
                    C->qual_scheme_size * sizeof(char));

            C->qual_scheme_lens =
                realloc_or_die(
                    C->qual_scheme_lens,
                    C->qual_scheme_size * sizeof(uint32_t));
        }

        last_base_qual = C->qual_scheme_vals[C->qual_scheme_count - 1] = min_qual;
        C->qual_scheme_lens[C->qual_scheme_count - 1] = C->chunk_len;
    }
    else {
        C->qual_scheme_lens[C->qual_scheme_count - 1] += C->chunk_len;
    }

    qualenc_set_base_qual(C->qualenc, last_base_qual);
}


static void quip_out_flush_chunk(quip_out_t* C)
{
    update_qual_scheme_guess(C);

    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    pthread_t id_thread, seq_thread, qual_thread;
    pthread_create_or_die(&id_thread,   &attr, id_compressor_thread,   (void*) C);
    pthread_create_or_die(&seq_thread,  &attr, seq_compressor_thread,  (void*) C);
    pthread_create_or_die(&qual_thread, &attr, qual_compressor_thread, (void*) C);

    size_t i;
    for (i = 0; i < C->chunk_len; ++i) {
        quip_out_add_readlen(C, C->chunk[i]->seq.n);
        C->id_bytes   += C->chunk[i]->id1.n;
        C->qual_bytes += C->chunk[i]->qual.n;
        C->seq_bytes  += C->chunk[i]->seq.n;
        C->buffered_bases += C->chunk[i]->seq.n;
        C->total_bases    += C->chunk[i]->seq.n;
    }

    C->buffered_reads += C->chunk_len;
    C->total_reads    += C->chunk_len;

    void* val;
    pthread_join_or_die(id_thread,   &val);
    pthread_join_or_die(seq_thread,  &val);
    pthread_join_or_die(qual_thread, &val);

    pthread_attr_destroy(&attr);

    C->chunk_len = 0;
}


int quip_out_readseq(quip_out_t* C, fastq_t* parser)
{
    if (C->buffered_bases > block_size) {
        quip_out_flush_block(C);
    }

    if (C->chunk_len == chunk_size) {
        quip_out_flush_chunk(C);
    }

    int r = fastq_next(parser, C->chunk[C->chunk_len]);

    if (r > 0) C->chunk_len++;
    return r;
}


void quip_out_addseq(quip_out_t* C, seq_t* seq)
{
    if (C->buffered_bases > block_size) {
        quip_out_flush_block(C);
    }

    if (C->chunk_len == chunk_size) {
        quip_out_flush_chunk(C);
    }

    fastq_copy_seq(C->chunk[C->chunk_len++], seq);
}


void quip_out_finish(quip_out_t* C)
{
    if (C->finished) return;
    if (C->chunk_len > 0) quip_out_flush_chunk(C);
    if (C->buffered_bases > 0) quip_out_flush_block(C);

    /* write an empty header to signify the end of the stream */
    write_uint32(C->writer, C->writer_data, 0);

    C->finished = true;
}


void quip_out_free(quip_out_t* C)
{
    if (!C->finished) quip_out_finish(C);

    size_t i;
    for (i = 0; i < chunk_size; ++i) {
        fastq_free_seq(C->chunk[i]);
    }

    idenc_free(C->idenc);
    qualenc_free(C->qualenc);
    assembler_free(C->assembler);
    free(C->readlen_vals);
    free(C->readlen_lens);
    free(C->qual_scheme_vals);
    free(C->qual_scheme_lens);
    free(C);
}


struct quip_in_t_
{
    /* sequence buffers */
    short_read_t* chunk[chunk_size];
    size_t chunk_len;
    size_t chunk_pos;

    /* function for writing compressed data */
    quip_reader_t reader;
    void* reader_data;

    /* algorithms to decompress ids, qualities, and sequences, resp. */
    idenc_t*        idenc;
    disassembler_t* disassembler;
    qualenc_t*      qualenc;

    /* compressed quality scores */
    uint8_t* qualbuf;
    size_t qualbuf_size, qualbuf_len, qualbuf_pos;

    /* compressed ids */
    uint8_t* idbuf;
    size_t idbuf_size, idbuf_len, idbuf_pos;

    /* compressed sequence */
    uint8_t* seqbuf;
    size_t seqbuf_size, seqbuf_len, seqbuf_pos;

    /* number of reads encoded in the buffers */
    uint32_t pending_reads;

    /* current block number */
    uint32_t block_num;

    /* expected block checksums */
    uint64_t exp_id_crc;
    uint64_t exp_seq_crc;
    uint64_t exp_qual_crc;

    /* observed block checksums */
    uint64_t id_crc;
    uint64_t seq_crc;
    uint64_t qual_crc;

    /* run length encoded read lengths */
    uint32_t* readlen_vals;
    uint32_t* readlen_lens;
    size_t readlen_count, readlen_size;

    size_t readlen_idx, readlen_off;

    /* run length encoded quality scheme guesses */
    char*     qual_scheme_vals;
    uint32_t* qual_scheme_lens;
    size_t qual_scheme_count, qual_scheme_size;
    size_t qual_scheme_idx, qual_scheme_off;

    bool end_of_stream;
};


static void* id_decompressor_thread(void* ctx)
{
    quip_in_t* D = (quip_in_t*) ctx;

    size_t cnt = D->pending_reads >= chunk_size ?
                    chunk_size : D->pending_reads;
    size_t i;
    for (i = 0; i < cnt; ++i) {
        idenc_decode(D->idenc, D->chunk[i]);
        D->id_crc = crc64_update(
            (uint8_t*) D->chunk[i]->id1.s,
            D->chunk[i]->id1.n, D->id_crc);
    }

    return NULL;
}


static void* seq_decompressor_thread(void* ctx)
{
    quip_in_t* D = (quip_in_t*) ctx;

    size_t readlen_idx = D->readlen_idx;
    size_t readlen_off = D->readlen_off;

    size_t n; /* read length */
    size_t cnt = D->pending_reads >= chunk_size ?
                    chunk_size : D->pending_reads;
    size_t i;

    for (i = 0; i < cnt; ) {
        n = D->readlen_vals[readlen_idx];
        if (++readlen_off >= D->readlen_lens[readlen_idx]) {
            readlen_off = 0;
            readlen_idx++;
        }

        disassembler_read(D->disassembler, D->chunk[i], n);
        D->seq_crc = crc64_update(
            (uint8_t*) D->chunk[i]->seq.s,
            D->chunk[i]->seq.n, D->seq_crc);
        ++i;
    }

    return NULL;
}


static void* qual_decompressor_thread(void* ctx)
{
    quip_in_t* D = (quip_in_t*) ctx;

    size_t readlen_idx = D->readlen_idx;
    size_t readlen_off = D->readlen_off;

    size_t qual_scheme_idx = D->qual_scheme_idx;
    size_t qual_scheme_off = D->qual_scheme_off;

    size_t n; /* read length */
    size_t cnt = D->pending_reads >= chunk_size ?
                    chunk_size : D->pending_reads;
    size_t i;
    for (i = 0; i < cnt; ++i) {
        n = D->readlen_vals[readlen_idx];
        if (++readlen_off >= D->readlen_lens[readlen_idx]) {
            readlen_off = 0;
            readlen_idx++;
        }

        if (++qual_scheme_off >= D->qual_scheme_lens[qual_scheme_idx]) {
            qual_scheme_off = 0;
            qual_scheme_idx++;
            if (qual_scheme_idx < D->qual_scheme_count) {
                qualenc_set_base_qual(D->qualenc, D->qual_scheme_vals[qual_scheme_idx]);
            }
        }

        qualenc_decode(D->qualenc, D->chunk[i], n);

        D->qual_crc = crc64_update(
            (uint8_t*) D->chunk[i]->qual.s,
            D->chunk[i]->qual.n, D->qual_crc);
    }

    return NULL;
}


static size_t qual_buf_reader(void* param, uint8_t* data, size_t size)
{
    quip_in_t* C = (quip_in_t*) param;

    size_t cnt = 0;
    while (cnt < size && C->qualbuf_pos < C->qualbuf_len) {
        data[cnt++] = C->qualbuf[C->qualbuf_pos++];
    }

    return cnt;
}

static size_t id_buf_reader(void* param, uint8_t* data, size_t size)
{
    quip_in_t* C = (quip_in_t*) param;

    size_t cnt = 0;
    while (cnt < size && C->idbuf_pos < C->idbuf_len) {
        data[cnt++] = C->idbuf[C->idbuf_pos++];
    }

    return cnt;
}


static size_t seq_buf_reader(void* param, uint8_t* data, size_t size)
{
    quip_in_t* C = (quip_in_t*) param;

    size_t cnt = 0;
    while (cnt < size && C->seqbuf_pos < C->seqbuf_len) {
        data[cnt++] = C->seqbuf[C->seqbuf_pos++];
    }

    return cnt;
}


quip_in_t* quip_in_alloc(quip_reader_t reader, void* reader_data)
{
    quip_in_t* D = malloc_or_die(sizeof(quip_in_t));

    D->reader = reader;
    D->reader_data = reader_data;

    size_t i;
    for (i = 0; i < chunk_size; ++i) {
        D->chunk[i] = fastq_alloc_seq();
    }
    D->chunk_len = 0;
    D->chunk_pos = 0;

    D->idenc   = idenc_alloc_decoder(id_buf_reader, (void*) D);
    D->disassembler = disassembler_alloc(seq_buf_reader, (void*) D);
    D->qualenc = qualenc_alloc_decoder(qual_buf_reader, (void*) D);

    D->qualbuf = NULL;
    D->qualbuf_size = 0;
    D->qualbuf_len  = 0;
    D->qualbuf_pos  = 0;

    D->idbuf = NULL;
    D->idbuf_size = 0;
    D->idbuf_len  = 0;
    D->idbuf_pos  = 0;

    D->seqbuf = NULL;
    D->seqbuf_size = 0;
    D->seqbuf_len  = 0;
    D->seqbuf_pos  = 0;

    D->pending_reads = 0;
    D->block_num = 0;

    D->readlen_size  = 1;
    D->readlen_count = 0;
    D->readlen_vals = malloc_or_die(D->readlen_size * sizeof(uint32_t));
    D->readlen_lens = malloc_or_die(D->readlen_size * sizeof(uint32_t));

    D->qual_scheme_size  = 1;
    D->qual_scheme_count = 0;
    D->qual_scheme_vals = malloc_or_die(D->qual_scheme_size * sizeof(char));
    D->qual_scheme_lens = malloc_or_die(D->qual_scheme_size * sizeof(uint32_t));

    D->id_crc   = D->exp_id_crc   = 0;
    D->seq_crc  = D->exp_seq_crc  = 0;
    D->qual_crc = D->exp_qual_crc = 0;

    D->end_of_stream = false;

    uint8_t header[7];
    if (D->reader(D->reader_data, header, 7) < 7 ||
        memcmp(quip_header_magic, header, 6) != 0) {
        fprintf(stderr, "Input is not a quip file.\n");
        exit(EXIT_FAILURE);
    }

    if (header[6] != quip_header_version) {
        fprintf(stderr, "Input is an old quip format --- an older version of quip is needed.\n");
        exit(EXIT_FAILURE);
    }

    return D;
}


void quip_in_free(quip_in_t* D)
{
    size_t i;
    for (i = 0; i < chunk_size; ++i) {
        fastq_free_seq(D->chunk[i]);
    }

    idenc_free(D->idenc);
    disassembler_free(D->disassembler);
    qualenc_free(D->qualenc);
    free(D->qualbuf);
    free(D->idbuf);
    free(D->seqbuf);
    free(D->readlen_vals);
    free(D->readlen_lens);
    free(D->qual_scheme_vals);
    free(D->qual_scheme_lens);
    free(D);
}


static void quip_in_read_block_header(quip_in_t* D)
{
    D->pending_reads = read_uint32(D->reader, D->reader_data);
    if (D->pending_reads == 0) {
        D->end_of_stream = true;
        return;
    }

    /* pending bases, which we don't need to know. */
    read_uint32(D->reader, D->reader_data);

    /* read run length encoded read lengths */
    uint32_t cnt = 0;
    uint32_t readlen_val, readlen_len;
    D->readlen_count = 0;
    while (cnt < D->pending_reads) {
        readlen_val = read_uint32(D->reader, D->reader_data);
        readlen_len = read_uint32(D->reader, D->reader_data);

        if (D->readlen_count >= D->readlen_size) {
            while (D->readlen_count >= D->readlen_size) D->readlen_size *= 2;
            D->readlen_vals = realloc_or_die(
                    D->readlen_vals, D->readlen_size * sizeof(uint32_t));
            D->readlen_lens = realloc_or_die(
                    D->readlen_lens, D->readlen_size * sizeof(uint32_t));
        }

        D->readlen_vals[D->readlen_count] = readlen_val;
        D->readlen_lens[D->readlen_count] = readlen_len;

        D->readlen_count++;
        cnt += readlen_len;
    }

    /* read run length encoded quality scheme guesses */
    cnt = 0;
    char     qual_scheme_val;
    uint32_t qual_scheme_len;
    D->qual_scheme_count = 0;
    while (cnt < D->pending_reads) {
        qual_scheme_val = (char) read_uint8(D->reader, D->reader_data);
        qual_scheme_len = read_uint32(D->reader, D->reader_data);

        if (D->qual_scheme_count >= D->qual_scheme_size) {
            while (D->qual_scheme_count >= D->qual_scheme_size) D->qual_scheme_size *= 2;
            D->qual_scheme_vals = realloc_or_die(
                    D->qual_scheme_vals, D->qual_scheme_size * sizeof(char));
            D->qual_scheme_lens = realloc_or_die(
                    D->qual_scheme_lens, D->qual_scheme_size * sizeof(uint32_t));
        }

        D->qual_scheme_vals[D->qual_scheme_count] = qual_scheme_val;
        D->qual_scheme_lens[D->qual_scheme_count] = qual_scheme_len;
        D->qual_scheme_count++;
        cnt += qual_scheme_len;
    }

    /* read id byte count */
    read_uint32(D->reader, D->reader_data); /* uncompressed bytes */
    uint32_t id_byte_cnt = read_uint32(D->reader, D->reader_data);
    if (id_byte_cnt > D->idbuf_size) {
        D->idbuf_size = id_byte_cnt;
        free(D->idbuf);
        D->idbuf = malloc_or_die(D->idbuf_size * sizeof(uint8_t));
    }
    D->exp_id_crc = read_uint64(D->reader, D->reader_data);

    /* read seq byte count */
    read_uint32(D->reader, D->reader_data); /* uncompressed bytes */
    uint32_t seq_byte_cnt = read_uint32(D->reader, D->reader_data);
    if (seq_byte_cnt > D->seqbuf_size) {
        D->seqbuf_size = seq_byte_cnt;
        free(D->seqbuf);
        D->seqbuf = malloc_or_die(D->seqbuf_size * sizeof(uint8_t));
    }
    D->exp_seq_crc = read_uint64(D->reader, D->reader_data);

    /* read qual byte count */
    read_uint32(D->reader, D->reader_data); /* uncompressed bytes */
    uint32_t qual_byte_cnt = read_uint32(D->reader, D->reader_data);
    if (qual_byte_cnt > D->qualbuf_size) {
        D->qualbuf_size = qual_byte_cnt;
        free(D->qualbuf);
        D->qualbuf = malloc_or_die(D->qualbuf_size * sizeof(uint8_t));
    }
    D->exp_qual_crc = read_uint64(D->reader, D->reader_data);

    if (D->pending_reads == 0 &&
             id_byte_cnt == 0 &&
            seq_byte_cnt == 0 &&
           qual_byte_cnt == 0)
    {
        D->end_of_stream = true;
        return;
    }

    /* read compressed data into buffers */
    D->idbuf_len = D->reader(D->reader_data, D->idbuf, id_byte_cnt);
    if (D->idbuf_len < id_byte_cnt) {
        fprintf(stderr, "Unexpected end of file.\n");
        exit(EXIT_FAILURE);
    }
    D->idbuf_pos = 0;

    D->seqbuf_len = D->reader(D->reader_data, D->seqbuf, seq_byte_cnt);
    if (D->seqbuf_len < seq_byte_cnt) {
        fprintf(stderr, "Unexpected end of file.\n");
        exit(EXIT_FAILURE);
    }
    D->seqbuf_pos = 0;

    D->qualbuf_len = D->reader(D->reader_data, D->qualbuf, qual_byte_cnt);
    if (D->qualbuf_len < qual_byte_cnt) {
        fprintf(stderr, "Unexpected end of file.\n");
        exit(EXIT_FAILURE);
    }
    D->qualbuf_pos = 0;

    D->readlen_idx = 0;
    D->readlen_off = 0;

    D->qual_scheme_off = 0;
    D->qual_scheme_idx = 0;
    /* skip over any quality schemes that were not used */
    while (D->qual_scheme_idx < D->qual_scheme_count - 1 &&
           D->qual_scheme_lens[D->qual_scheme_idx] == 0) {
        D->qual_scheme_idx++;
    }
    qualenc_set_base_qual(D->qualenc, D->qual_scheme_vals[D->qual_scheme_idx]);

    D->id_crc   = 0;
    D->seq_crc  = 0;
    D->qual_crc = 0;
    ++D->block_num;
}


seq_t* quip_in_read(quip_in_t* D)
{
    if (D->chunk_pos < D->chunk_len) {
        return D->chunk[D->chunk_pos++];
    }

    if (D->end_of_stream) return NULL;

    if (D->pending_reads == 0) {
        if (D->id_crc != D->exp_id_crc) {
            fprintf(stderr,
                "Warning: ID checksums in block %u do not match. "
                "ID data may be corrupt.\n", D->block_num);
        }

        if (D->seq_crc != D->exp_seq_crc) {
             fprintf(stderr,
                "Warning: Sequence checksums in block %u do not match. "
                "Sequence data may be corrupt.\n", D->block_num);
        }

        if (D->qual_crc != D->exp_qual_crc) {
              fprintf(stderr,
                "Warning: Quality checksums in block %u do not match. "
                "Quality data may be corrupt.\n", D->block_num);
        }

        quip_in_read_block_header(D);
        if (D->pending_reads == 0) return NULL;

        /* reset decoders */
        idenc_reset_decoder(D->idenc);
        idenc_start_decoder(D->idenc);

        disassembler_reset(D->disassembler);

        qualenc_reset_decoder(D->qualenc);
        qualenc_start_decoder(D->qualenc);
    }

    /* launch threads to decode a chunk of reads */
    pthread_t id_thread, seq_thread, qual_thread;

    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    pthread_create_or_die(&id_thread,   &attr, id_decompressor_thread,   (void*) D);
    pthread_create_or_die(&seq_thread,  &attr, seq_decompressor_thread,  (void*) D);
    pthread_create_or_die(&qual_thread, &attr, qual_decompressor_thread, (void*) D);

    pthread_join_or_die(id_thread,   NULL);
    pthread_join_or_die(seq_thread,  NULL);
    pthread_join_or_die(qual_thread, NULL);

    pthread_attr_destroy(&attr);

    D->chunk_len = D->pending_reads >= chunk_size ?
                    chunk_size : D->pending_reads;
    D->chunk_pos = 0;
 
    D->pending_reads -= D->chunk_len;
    size_t i;
    for (i = 0; i < D->chunk_len; ++i) {
        if (++D->readlen_off >= D->readlen_lens[D->readlen_idx]) {
            D->readlen_off = 0;
            D->readlen_idx++;
        }

        if (++D->qual_scheme_off >= D->qual_scheme_lens[D->qual_scheme_idx]) {
            D->qual_scheme_off = 0;
            D->qual_scheme_idx++;
        }
    }
    
    return D->chunk[D->chunk_pos++];
}


void quip_list(quip_reader_t reader, void* reader_data, quip_list_t* l)
{
    memset(l, 0, sizeof(quip_list_t));
    uint32_t block_reads;
    uint32_t readlen_count;
    uint32_t qual_scheme_count;
    uint64_t n;

    uint64_t block_bytes;

    uint8_t header[7];
    if (reader(reader_data, header, 7) < 7 ||
        memcmp(quip_header_magic, header, 6) != 0) {
        fprintf(stderr, "Input is not a quip file.\n");
        exit(EXIT_FAILURE);
    }

    if (header[6] != quip_header_version) {
        fprintf(stderr, "Input is an old quip format --- an older version of quip is needed.\n");
        exit(EXIT_FAILURE);
    }


    while (true) {
        /* read a block header */
        block_reads = read_uint32(reader, reader_data);
        l->header_bytes += 4;
        if (block_reads == 0) break;

        l->num_reads += block_reads;
        l->num_bases += read_uint32(reader, reader_data);
        l->num_blocks++;

        /* read lengths */
        readlen_count = 0;
        while (readlen_count < block_reads) {
            read_uint32(reader, reader_data); /* read length */
            readlen_count += read_uint32(reader, reader_data);
            l->header_bytes += 8;
        }

        qual_scheme_count = 0;
        while (qual_scheme_count < block_reads) {
            read_uint8(reader, reader_data); /* base_qual */
            qual_scheme_count += read_uint32(reader, reader_data);
            l->header_bytes += 6;
        }

        block_bytes = 0;

        /* id byte-count and checksum */
        l->id_bytes[0] += read_uint32(reader, reader_data);
        n = read_uint32(reader, reader_data);
        l->id_bytes[1] += n;
        block_bytes += n;
        read_uint64(reader, reader_data);

        /* sequence byte-count and checksum */
        l->seq_bytes[0] += read_uint32(reader, reader_data);
        n = read_uint32(reader, reader_data);
        l->seq_bytes[1] += n;
        block_bytes += n;
        read_uint64(reader, reader_data);

        /* quality scores byte-count and checksum */
        l->qual_bytes[0] += read_uint32(reader, reader_data);
        n = read_uint32(reader, reader_data);
        l->qual_bytes[1] += n;
        block_bytes += n;
        read_uint64(reader, reader_data);

        l->header_bytes += 52;

        /* seek past the compressed data */
        reader(reader_data, NULL, block_bytes);
    }
}

