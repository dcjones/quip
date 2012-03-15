/*
 * This file is part of quip.
 *
 * Copyright (c) 2012 by Daniel C. Jones <dcjones@cs.washington.edu>
 *
 */

#ifndef QUIP_QUIP
#define QUIP_QUIP

#include "quip_parse.h"
#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>

 /* All input and output in quip are performed by a simple callback interface.
  * The writer is responsible for writing compressed data, and has the takes
  * the following arguments:
  *     writer_data  :  user-supplied data used by the writer function
  *     data         :  data to be written
  *     size         :  number of bytes to write
  *
  * The reader function, used to read compressed data, works similarly but
  * with some more expectations. It's arguments:
  *     reader_data  :  user-supplied data used by the reader function
  *     data         :  buffer to write read data to
  *     size         :  number of bytes requested
  *
  * The reader function should read at most `size` bytes and write it
  * to the `data` buffer. It should then return the number of bytes
  * read. If the number returned is less than `size` it is taken to
  * mean than the end of the stream has been reached. 
  *
  * Furthermore, if data is NULL, no data need be written, but
  * the stream should advance by `size` bytes (or possible fewer, if
  * the end of the stream has been reached).
  *
  *
  */
typedef void   (*quip_writer_t) (void* writer_data, const uint8_t* data, size_t size);

typedef size_t (*quip_reader_t) (void* reader_data, uint8_t* data, size_t size);



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

/* Print a great deal of useless information while running. */
extern bool quip_verbose;

#endif


