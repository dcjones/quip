/*
 * This file is part of quip.
 *
 * Copyright (c) 2012 by Daniel C. Jones <dcjones@cs.washington.edu>
 *
 */


/*

Quip provides a unified interface for reading and writing short read
sequencing data is a variety of formats, no least of which, the quip format.

Currently the following formats are provided:

  * FASTQ
  * QUIP
  * SAM
  * BAM

Any format may be converted to any other, but conversions from SAM/BAM to
FASTQ will lose alignment information.

Additionally, our interpretation of the FASTQ format, excludes redundant
listing. If your FASTQ file looks like this:

  @unique_read_id_1234
  ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
  +unique_read_id_1234
  BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB

You will lose the second listing of 'unique_read_id_12234'. This entry will
be read as:

  @unique_read_id_1234
  ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
  +
  BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB

This is a feature, not a bug. We are doing you favor.

*/

#ifndef QUIP_QUIP
#define QUIP_QUIP

#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif


/* Simple encapsulated string. */
typedef struct
{
    unsigned char* s;    /* null-terminated string */
    size_t         n;    /* length of s */
    size_t         size; /* bytes allocated for s */
} str_t;

void str_init(str_t*);
void str_reserve(str_t*, size_t);
void str_free(str_t*);
void str_copy(str_t* dest, const str_t* src);
void str_copy_cstr(str_t* dest, const char* src, size_t n);


/* A series of CIGAR (edit) operations. */
typedef struct cigar_t_
{
    uint8_t*  ops;  /* cigar operations */
    uint32_t* lens; /* operation lengths */

    size_t n;    /* number of cigar ops */
    size_t size; /* size allocated */
} cigar_t;

void cigar_init(cigar_t*);
void cigar_reserve(cigar_t*, size_t);
void cigar_free(cigar_t*);


/* A read, either aligned or unaligned. */
typedef struct short_read_t_
{
    /* read ID */
    str_t id;

    /* Nucleotide sequence. A string over the alphabet [ACGTN] */
    str_t seq;

    /* Quality scores, of qual length as seq, from
     * ASCII chars in [33, 96] (that is, ['!', '`']) */
    str_t qual;

    /* SAM flags */
    uint32_t flags;

    /* Information for aligned reads. These must be set if
     * (flags & BAM_FUNMAP) == 0 (i.e., the read is not unmapped). */
    str_t    seqname;
    uint8_t  strand;
    uint32_t pos;
    cigar_t  cigar;

    /* auxiliary SAM fields */
    str_t aux;

} short_read_t;

void short_read_init(short_read_t*);
void short_read_free(short_read_t*);
void short_read_copy(short_read_t* dest, const short_read_t* src);


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


/*
 * Formats supported by quip.
 */
typedef enum quip_fmt_t_
{
    QUIP_FMT_UNDEFINED,
    QUIP_FMT_NULL,
    QUIP_FMT_FASTQ,
    QUIP_FMT_SAM,
    QUIP_FMT_BAM,
    QUIP_FMT_QUIP
} quip_fmt_t;


/*
 * Format specific options are specified using
 * combinations of these flags.
 */

/* Compress quip files without the de novo assembly step, which
 * is somewhat faster. */
#define QUIP_OPT_QUIP_ASSEMBLY_FREE 1

/* Output SAM files in BAM (compressed SAM) format. */
#define QUIP_OPT_SAM_BAM 1

typedef uint32_t quip_opt_t;


/*
 * Writing short read data in a variety of formats.
 */

typedef struct quip_out_t_ quip_out_t;

quip_out_t* quip_out_open(
              quip_writer_t writer,
              void*         writer_data,
              quip_fmt_t    format,
              quip_opt_t    opts);

void quip_out_close(quip_out_t*);

/*
 * Write one read to the output stream.
 */
void quip_write(quip_out_t*, short_read_t*);


/*
 * Reading short read data in a variety of formats.
 */

typedef struct quip_in_t_ quip_in_t;

quip_in_t* quip_in_open(
              quip_reader_t reader,
              void*         reader_data,
              quip_fmt_t    format,
              quip_opt_t    opts);


void quip_in_close(quip_in_t*);

/*
 * Get one read from the input stream. Output NULL if there are none.
 * The pointer returned is not guaranteed to be valid after
 * subsequent calls to quip_in_read, or to quip_in_close.
 */
short_read_t* quip_read(quip_in_t*);


/*
 * Read one read from the input stream and write it
 * to the output stream. Return false in the input
 * stream is empty, and true otherwise.
 */
bool quip_pipe(quip_in_t* in, quip_out_t* out);


/* Print a great deal of useless information while running. */
extern bool quip_verbose;


#ifdef __cplusplus
}
#endif
#endif


