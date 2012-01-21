
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


void quip_write_header(FILE* f)
{
    fwrite(quip_header_magic, 1, 6, f);
    fwrite(&quip_header_version, 1, 1, f);
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
    size_t* readlen_vals;
    size_t* readlen_lens;
    size_t readlen_count, readlen_size;

    /* general statistics */
    size_t total_reads;
    size_t total_bases;
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
    C->readlen_vals = malloc_or_die(C->readlen_size * sizeof(size_t));
    C->readlen_lens = malloc_or_die(C->readlen_size * sizeof(size_t));

    C->total_reads = 0;
    C->total_bases = 0;

    C->idenc     = idenc_alloc(id_buf_writer, (void*) C);
    C->qualenc   = qualenc_alloc(qual_buf_writer, (void*) C);
    C->assembler = assembler_alloc(seq_writer, (void*) C, assembler_k, aligner_k, quick);

    return C;
}

static void quip_comp_add_readlen(quip_compressor_t* C, size_t l)
{
    if (C->readlen_count == 0 || C->readlen_vals[C->readlen_count - 1] != l) {
        if (C->readlen_count >= C->readlen_size) {
            while (C->readlen_count >= C->readlen_size) {
                C->readlen_size *= 2;
            }

            C->readlen_vals = realloc_or_die(C->readlen_vals, C->readlen_size * sizeof(size_t));
            C->readlen_lens = realloc_or_die(C->readlen_vals, C->readlen_size * sizeof(size_t));
        }

        C->readlen_vals[C->readlen_count] = l;
        C->readlen_lens[C->readlen_count] = 1;
        C->readlen_count++;
    }
    else C->readlen_lens[C->readlen_count - 1]++;
}


void quip_comp_addseq(quip_compressor_t* C, seq_t* seq)
{
    if (C->buffered_bases > block_size) {
        quip_comp_flush(C);
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

    /* write ids */
    idenc_flush(C->idenc);
    C->writer(C->writer_data, C->idbuf, C->idbuf_len);
    if (verbose) {
        fprintf(stderr, "\tid: %zu / %zu (%0.2f%%)\n",
                C->idbuf_len, C->id_bytes,
                100.0 * (double) C->idbuf_len / (double) C->id_bytes);
    }
    C->idbuf_len = 0;

    /* write sequences */
    assembler_assemble(C->assembler);
    C->writer(C->writer_data, C->seqbuf, C->seqbuf_len);
    if (verbose) {
        fprintf(stderr, "\tseq: %zu / %zu (%0.2f%%)\n",
                C->seqbuf_len, C->seq_bytes,
                100.0 * (double) C->seqbuf_len / (double) C->seq_bytes);
    }
    C->seqbuf_len = 0;

    /* write qualities */
    qualenc_flush(C->qualenc);
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


void quip_comp_free(quip_compressor_t* C)
{
    quip_comp_flush(C);
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


