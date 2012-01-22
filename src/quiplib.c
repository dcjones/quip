
#include "quip.h"
#include "assembler.h"
#include "misc.h"
#include "qualenc.h"
#include "idenc.h"
#include <stdint.h>
#include <string.h>

const uint8_t quip_header_magic[6] = 
    {0xff, 'Q', 'U', 'I', 'P', 0x00};

const uint8_t quip_header_version = 0x01;

const size_t assembler_k = 25;
const size_t aligner_k = 15;

/* approximate number of bases per block */
const size_t block_size = 100000000;

bool verbose = true;


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


struct quip_compressor_t_
{
    /* function for writing compressed data */
    quip_writer_t writer;
    void* writer_data;

    /* algorithms to compress ids, qualities, and sequences, resp. */
    idenc_t*     idenc;
    qualenc_t*   qualenc;
    assembler_t* assembler;

    /* compressed quality scores */
    uint8_t* qualbuf;
    size_t qualbuf_size, qualbuf_len;

    /* compressed ids */
    uint8_t* idbuf;
    size_t idbuf_size, idbuf_len;

    /* compressed sequence */
    uint8_t* seqbuf;
    size_t seqbuf_size, seqbuf_len;

    /* number of bases encoded in the current block */
    size_t buffered_reads;
    size_t buffered_bases;

    /* for compression statistics */
    size_t id_bytes;
    size_t qual_bytes;
    size_t seq_bytes;

    /* run length encoded read lengths */
    uint32_t* readlen_vals;
    uint32_t* readlen_lens;
    size_t readlen_count, readlen_size;

    /* general statistics */
    uint64_t total_reads;
    uint64_t total_bases;
};


static void qual_buf_writer(void* param, const uint8_t* data, size_t size)
{
    quip_compressor_t* C = (quip_compressor_t*) param;

    if (C->qualbuf_size < C->qualbuf_len + size) {
        while (C->qualbuf_size < C->qualbuf_len + size) {
            C->qualbuf_size *= 2;
        }

        C->qualbuf = realloc_or_die(C->qualbuf, C->qualbuf_size * sizeof(uint8_t));
    }

    memcpy(C->qualbuf + C->qualbuf_len, data, size);
    C->qualbuf_len += size;
}


static void id_buf_writer(void* param, const uint8_t* data, size_t size)
{
    quip_compressor_t* C = (quip_compressor_t*) param;

    if (C->idbuf_size < C->idbuf_len + size) {
        while (C->idbuf_size < C->idbuf_len + size) {
            C->idbuf_size *= 2;
        }

        C->idbuf = realloc_or_die(C->idbuf, C->idbuf_size * sizeof(uint8_t));
    }

    memcpy(C->idbuf + C->idbuf_len, data, size);
    C->idbuf_len += size;
}


static void seq_writer(void* param, const uint8_t* data, size_t size)
{
    quip_compressor_t* C = (quip_compressor_t*) param;

    if (C->seqbuf_size < C->seqbuf_len + size) {
        while (C->seqbuf_size < C->seqbuf_len + size) {
            C->seqbuf_size *= 2;
        }

        C->seqbuf = realloc_or_die(C->seqbuf, C->seqbuf_size * sizeof(uint8_t));
    }

    memcpy(C->seqbuf + C->seqbuf_len, data, size);
    C->seqbuf_len += size;
}



quip_compressor_t* quip_comp_alloc(quip_writer_t writer, void* writer_data, bool quick)
{
    quip_compressor_t* C = malloc_or_die(sizeof(quip_compressor_t));
    C->writer = writer;
    C->writer_data = writer_data;

    C->qualbuf_size = 1024;
    C->qualbuf_len = 0;
    C->qualbuf = malloc_or_die(C->qualbuf_size * sizeof(uint8_t));

    C->idbuf_size = 1024;
    C->idbuf_len = 0;
    C->idbuf = malloc_or_die(C->idbuf_size * sizeof(uint8_t));

    C->seqbuf_size = 1024;
    C->seqbuf_len = 0;
    C->seqbuf = malloc_or_die(C->seqbuf_size * sizeof(uint8_t));

    C->buffered_reads = 0;
    C->buffered_bases = 0;

    C->id_bytes   = 0;
    C->qual_bytes = 0;
    C->seq_bytes  = 0;

    C->readlen_size  = 1;
    C->readlen_count = 0;
    C->readlen_vals = malloc_or_die(C->readlen_size * sizeof(uint32_t));
    C->readlen_lens = malloc_or_die(C->readlen_size * sizeof(uint32_t));

    C->total_reads = 0;
    C->total_bases = 0;

    C->idenc     = idenc_alloc_encoder(id_buf_writer, (void*) C);
    C->qualenc   = qualenc_alloc_encoder(qual_buf_writer, (void*) C);
    C->assembler = assembler_alloc(seq_writer, (void*) C, assembler_k, aligner_k, quick);

    /* write header */
    C->writer(C->writer_data, quip_header_magic, 6);
    C->writer(C->writer_data, &quip_header_version, 1);


    return C;
}

static void quip_comp_add_readlen(quip_compressor_t* C, size_t l)
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


void quip_comp_flush_block(quip_compressor_t* C)
{
    if (verbose) {
        fprintf(stderr, "writing a block of %zu compressed bases...\n", C->buffered_bases);
    }

    /* write the number of encoded reads in the block */
    write_uint32(C->writer, C->writer_data, C->buffered_reads);

    /* write the (run length encoded) read lengths. */
    size_t i;
    for (i = 0; i < C->readlen_count; ++i) {
        write_uint32(C->writer, C->writer_data, C->readlen_vals[i]);
        write_uint32(C->writer, C->writer_data, C->readlen_lens[i]);
    }

    /* finish coding */
    idenc_flush(C->idenc);
    assembler_assemble(C->assembler);
    qualenc_flush(C->qualenc);

    /* write the size of each chunk about to be written */
    write_uint32(C->writer, C->writer_data, C->idbuf_len);
    write_uint32(C->writer, C->writer_data, C->seqbuf_len);
    write_uint32(C->writer, C->writer_data, C->qualbuf_len);

    /* write compressed ids */
    C->writer(C->writer_data, C->idbuf, C->idbuf_len);
    if (verbose) {
        fprintf(stderr, "\tid: %zu / %zu (%0.2f%%)\n",
                C->idbuf_len, C->id_bytes,
                100.0 * (double) C->idbuf_len / (double) C->id_bytes);
    }
    C->idbuf_len = 0;

    /* write compressed sequences */
    C->writer(C->writer_data, C->seqbuf, C->seqbuf_len);
    if (verbose) {
        fprintf(stderr, "\tseq: %zu / %zu (%0.2f%%)\n",
                C->seqbuf_len, C->seq_bytes,
                100.0 * (double) C->seqbuf_len / (double) C->seq_bytes);
    }
    C->seqbuf_len = 0;

    /* write compressed qualities */
    C->writer(C->writer_data, C->qualbuf, C->qualbuf_len);
    if (verbose) {
        fprintf(stderr, "\tqual: %zu / %zu (%0.2f%%)\n",
                C->qualbuf_len, C->qual_bytes,
                100.0 * (double) C->qualbuf_len / (double) C->qual_bytes);
    }
    C->qualbuf_len = 0;

    /* reset */
    C->buffered_reads = 0;
    C->buffered_bases = 0;
    C->id_bytes   = 0;
    C->qual_bytes = 0;
    C->seq_bytes  = 0;
    C->readlen_count = 0;
}


void quip_comp_addseq(quip_compressor_t* C, seq_t* seq)
{
    if (C->buffered_bases > block_size) {
        quip_comp_flush_block(C);
    }

    idenc_encode(C->idenc, seq);
    qualenc_encode(C->qualenc, seq);
    assembler_add_seq(C->assembler, seq->seq.s, seq->seq.n);

    C->buffered_reads++;
    C->buffered_bases += seq->seq.n;

    C->id_bytes   += seq->id1.n;
    C->qual_bytes += seq->qual.n;
    C->seq_bytes  += seq->seq.n;

    quip_comp_add_readlen(C, seq->seq.n);
    C->total_reads++;
    C->total_bases += seq->qual.n;
}


void quip_comp_flush(quip_compressor_t* C)
{
    if (C->buffered_bases > 0) quip_comp_flush_block(C);

    /* write an empty header to signify the end of the stream */
    write_uint32(C->writer, C->writer_data, 0);
    write_uint32(C->writer, C->writer_data, 0);
    write_uint32(C->writer, C->writer_data, 0);
    write_uint32(C->writer, C->writer_data, 0);

    /* write the total number of bases */
    write_uint64(C->writer, C->writer_data, C->total_bases);

    /* write the total number of reads */
    write_uint64(C->writer, C->writer_data, C->total_reads);
}



void quip_comp_free(quip_compressor_t* C)
{
    idenc_free(C->idenc);
    qualenc_free(C->qualenc);
    assembler_free(C->assembler);
    free(C->qualbuf);
    free(C->idbuf);
    free(C->seqbuf);
    free(C->readlen_vals);
    free(C->readlen_lens);
    free(C);
}



struct quip_decompressor_t_
{
    /* function for writing compressed data */
    quip_reader_t reader;
    void* reader_data;

    /* algorithms to decompress ids, qualities, and sequences, resp. */
    idenc_t*     idenc;
    qualenc_t*   qualenc;
    // TODO: we need a special dissasembler */
    /*assembler_t* assembler;*/

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
    size_t pending_reads;

    /* run length encoded read lengths */
    uint32_t* readlen_vals;
    uint32_t* readlen_lens;
    size_t readlen_count, readlen_size;

    size_t readlen_idx, readlen_off;
};


static size_t qual_buf_reader(void* param, uint8_t* data, size_t size)
{
    quip_decompressor_t* C = (quip_decompressor_t*) param;

    size_t cnt = 0;
    while (cnt < size && C->qualbuf_pos < C->qualbuf_len) {
        data[cnt++] = C->qualbuf[C->qualbuf_pos++];
    }

    return cnt;
}

static size_t id_buf_reader(void* param, uint8_t* data, size_t size)
{
    quip_decompressor_t* C = (quip_decompressor_t*) param;

    size_t cnt = 0;
    while (cnt < size && C->idbuf_pos < C->idbuf_len) {
        data[cnt++] = C->idbuf[C->idbuf_pos++];
    }

    return cnt;
}


static size_t seq_buf_reader(void* param, uint8_t* data, size_t size)
{
    quip_decompressor_t* C = (quip_decompressor_t*) param;

    size_t cnt = 0;
    while (cnt < size && C->seqbuf_pos < C->seqbuf_len) {
        data[cnt++] = C->seqbuf[C->seqbuf_pos++];
    }

    return cnt;
}


quip_decompressor_t* quip_decomp_alloc(quip_reader_t reader, void* reader_data)
{
    quip_decompressor_t* D = malloc_or_die(sizeof(quip_decompressor_t));

    D->reader = reader;
    D->reader_data = reader_data;

    D->qualenc = qualenc_alloc_decoder(qual_buf_reader, (void*) D);
    D->idenc   = idenc_alloc_decoder(id_buf_reader, (void*) D);

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

    D->readlen_size  = 1;
    D->readlen_count = 0;
    D->readlen_vals = malloc_or_die(D->readlen_size * sizeof(uint32_t));
    D->readlen_lens = malloc_or_die(D->readlen_size * sizeof(uint32_t));

    uint8_t header[7];
    D->reader(D->reader_data, header, 7);
    if (memcmp(quip_header_magic, header, 6) != 0) {
        fprintf(stderr, "Input is not a quip file.\n");
        exit(EXIT_FAILURE);
    }

    if (header[6] != quip_header_version) {
        fprintf(stderr, "Input is an old quip format --- an older version of quip is needed.\n");
        exit(EXIT_FAILURE);
    }

    return D;
}


void quip_decomp_free(quip_decompressor_t* D)
{
    idenc_free(D->idenc);
    qualenc_free(D->qualenc);
    free(D->qualbuf);
    free(D->idbuf);
    free(D->seqbuf);
    free(D);
}


static void quip_decomp_read_block_header(quip_decompressor_t* D)
{
    D->pending_reads = read_uint32(D->reader, D->reader_data);

    /* read run length encoded read lengths */
    uint32_t cnt = 0;
    uint32_t val, len;
    D->readlen_count = 0;
    while (cnt < D->pending_reads) {
        val = read_uint32(D->reader, D->reader_data);
        len = read_uint32(D->reader, D->reader_data);

        if (D->readlen_count >= D->readlen_size) {
            while (D->readlen_count >= D->readlen_size) D->readlen_size *= 2;
            D->readlen_vals = realloc_or_die(
                    D->readlen_vals, D->readlen_size * sizeof(uint32_t));
            D->readlen_lens = realloc_or_die(
                    D->readlen_lens, D->readlen_size * sizeof(uint32_t));
        }

        D->readlen_vals[D->readlen_count] = val;
        D->readlen_lens[D->readlen_count] = len;

        D->readlen_count++;
        cnt += len;
    }

    /* read id byte count */
    uint32_t id_byte_cnt = read_uint32(D->reader, D->reader_data);
    if (id_byte_cnt > D->idbuf_size) {
        D->idbuf_size = id_byte_cnt;
        free(D->idbuf);
        D->idbuf = malloc_or_die(D->idbuf_size * sizeof(uint8_t));
    }

    /* read seq byte count */
    uint32_t seq_byte_cnt = read_uint32(D->reader, D->reader_data);
    if (seq_byte_cnt > D->seqbuf_size) {
        D->seqbuf_size = seq_byte_cnt;
        free(D->seqbuf);
        D->seqbuf = malloc_or_die(D->seqbuf_size * sizeof(uint8_t));
    }

    /* read qual byte count */
    uint32_t qual_byte_cnt = read_uint32(D->reader, D->reader_data);
    if (qual_byte_cnt > D->qualbuf_size) {
        D->qualbuf_size = qual_byte_cnt;
        free(D->qualbuf);
        D->qualbuf = malloc_or_die(D->qualbuf_size * sizeof(uint8_t));
    }

    /* read compressed data into buffers */
    D->idbuf_len = D->reader(D->reader_data, D->idbuf, id_byte_cnt);
    if (D->idbuf_len < id_byte_cnt) {
        fprintf(stderr, "Unexpected end of file.\n");
        exit(EXIT_FAILURE);
    }

    D->seqbuf_len = D->reader(D->reader_data, D->seqbuf, seq_byte_cnt);
    if (D->seqbuf_len < seq_byte_cnt) {
        fprintf(stderr, "Unexpected end of file.\n");
        exit(EXIT_FAILURE);
    }

    D->qualbuf_len = D->reader(D->reader_data, D->qualbuf, qual_byte_cnt);
    if (D->qualbuf_len < qual_byte_cnt) {
        fprintf(stderr, "Unexpected end of file.\n");
        exit(EXIT_FAILURE);
    }

    D->readlen_idx = 0;
    D->readlen_off = 0;
}


bool quip_decomp_read(quip_decompressor_t* D, seq_t* seq)
{
    if (D->pending_reads == 0) {
        quip_decomp_read_block_header(D);
        if (D->pending_reads == 0) return false;
    }

    /* read length */
    size_t n = D->readlen_vals[D->readlen_idx];
    if (++D->readlen_off >= D->readlen_lens[D->readlen_idx]) {
        D->readlen_off = 0;
        D->readlen_idx++;
    }

    idenc_decode(D->idenc, seq);

    // TODO: decode sequence

    qualenc_decode(D->qualenc, seq, n);

    D->pending_reads--;

    return true;
}


