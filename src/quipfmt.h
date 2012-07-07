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
                    quip_writer_t     writer,
                    void*             writer_data,
                    quip_opt_t        opts,
                    const quip_aux_t* aux,
                    const seqmap_t*   ref);

void quip_quip_out_close(quip_quip_out_t*);
void quip_quip_write(quip_quip_out_t*, short_read_t*);


typedef struct quip_quip_in_t_ quip_quip_in_t;

quip_quip_in_t* quip_quip_in_open(
                    quip_reader_t   reader,
                    void*           reader_data,
                    quip_opt_t      opts,
                    const seqmap_t* ref);

void quip_quip_in_close(quip_quip_in_t*);
void quip_quip_get_aux(quip_quip_in_t*, quip_aux_t*);
short_read_t* quip_quip_read(quip_quip_in_t*);


/* Efficiently determine the number of reads summary information
   about a compressed stream. */
typedef struct quip_list_t_
{
    uint64_t num_blocks;
    uint64_t num_bases;
    uint64_t num_reads;

    /* leading extra bytes (e.g. SAM header) */
    quip_fmt_t lead_fmt;
    uint64_t   lead_bytes;

    /* the uncompressed (0) and compressed (1) byte counts */
    uint64_t id_bytes[2];
    uint64_t aux_bytes[2];
    uint64_t seq_bytes[2];
    uint64_t qual_bytes[2];
    uint64_t header_bytes;

} quip_list_t;

void quip_list(quip_reader_t, void* reader_data, quip_list_t*);
void quip_list_file(FILE* file, quip_list_t*);

#endif

