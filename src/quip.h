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
    uint32_t       n;    /* length of s */
    uint32_t       size; /* bytes allocated for s */
} str_t;

void str_init(str_t*);
void str_reserve(str_t*, size_t);
void str_reserve_extra(str_t*, size_t);
void str_append(str_t*, const str_t*);
void str_append_cstr(str_t*, const char*);
void str_append_char(str_t*, char);
void str_free(str_t*);
void str_copy(str_t* dest, const str_t* src);
void str_copy_cstr(str_t* dest, const char* src, size_t n);
void str_memcpy(str_t* dest, const uint8_t* src, size_t n);


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
void cigar_copy(cigar_t*, const cigar_t*);


/* Sam optional fields. */
typedef struct samopt_t_
{
    unsigned char key[2];
    unsigned char type;
    str_t*        data;
} samopt_t;

typedef struct samopt_table_t_ samopt_table_t;
samopt_table_t* samopt_table_alloc();
void            samopt_table_free(samopt_table_t*);
void            samopt_table_clear(samopt_table_t*);
samopt_t*       samopt_table_get(samopt_table_t* M, const unsigned char key[2]);
size_t          samopt_table_size(const samopt_table_t* M);
void            samopt_table_copy(samopt_table_t* dest, const samopt_table_t* src);
size_t          samopt_table_bytes(const samopt_table_t*);
uint64_t        samopt_table_crc64_update(const samopt_table_t*, uint64_t crc);
void            samopt_table_dump_sorted(const samopt_table_t*, samopt_t** opts);

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
    uint8_t  map_qual;
    cigar_t  cigar;

    str_t    mate_seqname;
    uint32_t mate_pos;
    int32_t  tlen;

    /* optional SAM fields */
    samopt_table_t* aux;

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



void writer(void* param, const uint8_t* data, size_t datalen);
size_t quip_file_reader(void* param, uint8_t* data, size_t datalen);


/* Some hand read/write functions. */
void write_uint8(quip_writer_t writer, void* writer_data, uint8_t x);
void write_uint32(quip_writer_t writer, void* writer_data, uint32_t x);
void write_uint64(quip_writer_t writer, void* writer_data, uint64_t x);
uint8_t read_uint8(quip_reader_t reader, void* reader_data);
uint32_t read_uint32(quip_reader_t reader, void* reader_data);
uint64_t read_uint64(quip_reader_t reader, void* reader_data);


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

/* Compress quip files with a de novo assembly step: better compression
 * at the cost of compression and decompression speed. */
#define QUIP_OPT_QUIP_ASSEMBLY 1

/* Output SAM files in BAM (compressed SAM) format. */
#define QUIP_OPT_SAM_BAM 1

typedef uint32_t quip_opt_t;

typedef struct seqmap_t_ seqmap_t;

/* 
 * Certain formats carry with the auxiliary data.
 * This is a container to facilitate passing this
 * around.
 */
typedef struct quip_aux_t_
{
    quip_fmt_t fmt;
    str_t data;
} quip_aux_t;


/*
 * Writing short read data in a variety of formats.
 */

typedef struct quip_out_t_ quip_out_t;

quip_out_t* quip_out_open(
              quip_writer_t     writer,
              void*             writer_data,
              quip_fmt_t        format,
              quip_opt_t        opts,
              const quip_aux_t* aux,
              const seqmap_t*   ref);


quip_out_t* quip_out_open_file(
              FILE*             file,
              quip_fmt_t        format,
              quip_opt_t        opts,
              const quip_aux_t* aux,
              const seqmap_t*   ref);


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
              quip_reader_t   reader,
              void*           reader_data,
              quip_fmt_t      format,
              quip_opt_t      opts,
              const seqmap_t* ref);

quip_in_t* quip_in_open_file(
              FILE*           file,
              quip_fmt_t      format,
              quip_opt_t      opts,
              const seqmap_t* ref);

void quip_in_close(quip_in_t*);

void quip_get_aux(quip_in_t*, quip_aux_t*);

/*
 * Get one read from the input stream. Output NULL if there are none.
 * The pointer returned should not be freed by the caller, and is not
 * guaranteed to be valid after subsequent calls to quip_in_read, or
 * to quip_in_close.
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

/* Program name. */
extern const char* quip_prog_name;

/* Input file name. */
extern const char* quip_in_fname;

/* Number of reads used for assembly. */
extern size_t quip_assembly_n;

/* Exit with cleanup. */
void quip_abort();

/* Fatal errors. */
void quip_error(const char* fmt, ...);

/* Nonfatal errors. */
void quip_warning(const char* fmt, ...);



#ifdef __cplusplus
}
#endif
#endif


