/*
 * This file is part of quip.
 *
 * Copyright (c) 2012 by Daniel C. Jones <dcjones@cs.washington.edu>
 *
 */

#ifndef QUIP_QUIPFMT
#define QUIP_QUIPFMT

#include "quip.h"


typedef struct quip_quip_out_t_ quip_quip_out_t;

quip_quip_out_t* quip_quip_out_open(
                    quip_writer_t writer,
                    void*         writer_data,
                    quip_opt_t    opts);

void quip_quip_out_close(quip_quip_out_t*);
void quip_quip_write(quip_quip_out_t*, short_read_t*);


typedef struct quip_quip_in_t_ quip_quip_in_t;

quip_quip_in_t* quip_quip_in_open(
                    quip_reader_t reader,
                    void*         reader_data,
                    quip_opt_t opts);

void quip_quip_in_close(quip_quip_in_t*);
short_read_t* quip_quip_read(quip_quip_in_t*);


#if 0


/* Compression */
typedef struct quip_out_t_ quip_out_t;


/* Create a new compressor with the given stream of raw fastq data. */
quip_out_t* quip_out_alloc(quip_writer_t, void* writer_data, bool quick);

/* Deallocate a compressor, once finished with it.
 * This will flush any buffered data to the output stream. */
void quip_out_free(quip_out_t*);

/* Compress a single read. */
void quip_out_addseq(quip_out_t*, seq_t*);

/* Compress a single read, read from a parser (slightly more efficient than addseq) */
int quip_out_readseq(quip_out_t*, fastq_t*);

/* Finish compressing reads, flushing all buffered data,
 * and writing a end of stream marker. */
void quip_out_finish(quip_out_t*);



/* Decompression */
typedef struct quip_in_t_ quip_in_t;

/* Create a new decompressor with the given stream of compressed data in the quip format. */
quip_in_t* quip_in_alloc(quip_reader_t, void* reader_data);

/* Deallocate a decompressor, once finished with it. */
void quip_in_free(quip_in_t*);

/* Decompress one read from the stream. The pointer returned
 * should not be freed and is not guaranteed to be valid after
 * subsequent calls to quip_in_read, or calls to quip_in_free.
 *
 * If there are no more reads in the stream, NULL is returned.
 */
seq_t* quip_in_read(quip_in_t*);


/* Efficiently determine the number of reads summary information
   about a compressed stream. */
typedef struct quip_list_t_
{
    uint64_t num_blocks;
    uint64_t num_bases;
    uint64_t num_reads;

    /* the uncompressed (0) and compressed (1) byte counts */
    uint64_t id_bytes[2];
    uint64_t seq_bytes[2];
    uint64_t qual_bytes[2];
    uint64_t header_bytes;

} quip_list_t;

void quip_list(quip_reader_t, void* reader_data, quip_list_t*);

/* Test the integrity of a stream. */
void quip_test(quip_reader_t, void* reader_data); 

#endif

#endif

