
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


void quip_write_header(FILE* f)
{
    fwrite(quip_header_magic, 1, 6, f);
    fwrite(&quip_header_version, 1, 1, f);
}


struct quip_compressor_t_
{
    /* function for writing compressed data */
    quip_block_writer_t writer;
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


quip_compressor_t* quip_comp_alloc(quip_block_writer_t writer, void* writer_data)
{
    quip_compressor_t* C = malloc_or_die(sizeof(quip_compressor_t));
    C->writer = writer;
    C->writer_data = writer_data;

    C->idenc     = idenc_alloc(id_buf_writer, (void*) C);
    C->qualenc   = qualenc_alloc(qual_buf_writer, (void*) C);
    C->assembler = assembler_alloc(assembler_k, aligner_k);

    C->qualbuf_size = 1024;
    C->qualbuf_len = 0;
    C->qualbuf = malloc_or_die(C->qualbuf_size * sizeof(uint8_t));

    C->idbuf_size = 1024;
    C->idbuf_len = 0;
    C->idbuf = malloc_or_die(C->idbuf_size * sizeof(uint8_t));

    C->buffered_bases = 0;

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
}


void quip_comp_flush(quip_compressor_t* C)
{
    idenc_flush(C->idenc);
    C->writer(C->writer_data, C->idbuf, C->idbuf_len);
    C->idbuf_len = 0;

    qualenc_flush(C->qualenc);
    C->writer(C->writer_data, C->qualbuf, C->qualbuf_len);
    C->qualbuf_len = 0;

    assembler_assemble(C->assembler, C->writer, C->writer_data);
    assembler_flush(C->assembler);

    C->buffered_bases = 0;
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



