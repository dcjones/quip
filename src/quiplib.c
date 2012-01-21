
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

    /* algorithms to compress ids, qualities, and sequences, reps. */
    idenc_t*     idenc;
    qualenc_t*   qualenc;
    assembler_t* assembler;

    /* compressed quality scores */
    uint8_t* qualbuf;
    size_t qualbuf_size, qualbuf_len;

    /* compressed ids */
    uint8_t* idbuf;
    size_t idbuf_size, idbuf_len;

    /* number of bases encoded in the current block */
    size_t buffered_bases;

    /* for compression statistics */
    size_t id_bytes;
    size_t qual_bytes;
    size_t seq_bytes, seq_comp_bytes;
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
    C->seq_comp_bytes += size;
    C->writer(C->writer_data, data, size);
}



quip_compressor_t* quip_comp_alloc(quip_writer_t writer, void* writer_data, bool quick)
{
    quip_compressor_t* C = malloc_or_die(sizeof(quip_compressor_t));
    C->writer = writer;
    C->writer_data = writer_data;

    C->idenc     = idenc_alloc(id_buf_writer, (void*) C);
    C->qualenc   = qualenc_alloc(qual_buf_writer, (void*) C);
    C->assembler = assembler_alloc(seq_writer, (void*) C, assembler_k, aligner_k, quick);

    C->qualbuf_size = 1024;
    C->qualbuf_len = 0;
    C->qualbuf = malloc_or_die(C->qualbuf_size * sizeof(uint8_t));

    C->idbuf_size = 1024;
    C->idbuf_len = 0;
    C->idbuf = malloc_or_die(C->idbuf_size * sizeof(uint8_t));

    C->buffered_bases = 0;

    C->id_bytes   = 0;
    C->qual_bytes = 0;
    C->seq_bytes  = 0;

    return C;
}


void quip_comp_addseq(quip_compressor_t* C, seq_t* seq)
{
    if (C->buffered_bases > block_size) {
        quip_comp_flush(C);
    }

    idenc_encode(C->idenc, seq);
    qualenc_encode(C->qualenc, seq);
    assembler_add_seq(C->assembler, seq->seq.s, seq->seq.n);

    C->buffered_bases += seq->seq.n;

    C->id_bytes   += seq->id1.n;
    C->qual_bytes += seq->qual.n;
    C->seq_bytes  += seq->seq.n;
}


void quip_comp_flush(quip_compressor_t* C)
{
    if (verbose) {
        fprintf(stderr, "writing a block of %zu compressed bases...\n", C->buffered_bases);
    }

    assembler_assemble(C->assembler);
    if (verbose) {
        fprintf(stderr, "\tseq: %zu / %zu (%0.2f%%)\n",
                C->seq_comp_bytes, C->seq_bytes,
                100.0 * (double) C->seq_comp_bytes / (double) C->seq_bytes);
    }
    C->seq_comp_bytes = 0;


    idenc_flush(C->idenc);
    C->writer(C->writer_data, C->idbuf, C->idbuf_len);
    if (verbose) {
        fprintf(stderr, "\tid: %zu / %zu (%0.2f%%)\n",
                C->idbuf_len, C->id_bytes,
                100.0 * (double) C->idbuf_len / (double) C->id_bytes);
    }
    C->idbuf_len = 0;

    qualenc_flush(C->qualenc);
    C->writer(C->writer_data, C->qualbuf, C->qualbuf_len);
    if (verbose) {
        fprintf(stderr, "\tqual: %zu / %zu (%0.2f%%)\n",
                C->qualbuf_len, C->qual_bytes,
                100.0 * (double) C->qualbuf_len / (double) C->qual_bytes);
    }
    C->qualbuf_len = 0;


    C->buffered_bases = 0;
    C->id_bytes   = 0;
    C->qual_bytes = 0;
    C->seq_bytes  = 0;
}


void quip_comp_free(quip_compressor_t* C)
{
    quip_comp_flush(C);
    idenc_free(C->idenc);
    qualenc_free(C->qualenc);
    assembler_free(C->assembler);
    free(C->qualbuf);
    free(C->idbuf);
    free(C);
}



