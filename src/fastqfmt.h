/*
 * This file is part of quip.
 *
 * Copyright (c) 2012 by Daniel C. Jones <dcjones@cs.washington.edu>
 *
 */

#ifndef QUIP_FASTQFMT
#define QUIP_FASTQFMT

#include "quip.h"

typedef struct quip_fastq_out_t_ quip_fastq_out_t;

quip_fastq_out_t* quip_fastq_out_open(
                    quip_writer_t writer,
                    void*         writer_data,
                    quip_opt_t    opts);

void quip_fastq_out_close(quip_fastq_out_t*);
void quip_fastq_write(quip_fastq_out_t*, short_read_t*);


typedef struct quip_fastq_in_t_ quip_fastq_in_t;

quip_fastq_in_t* quip_fastq_in_open(
                    quip_reader_t reader,
                    void*         reader_data,
                    quip_opt_t opts);

void quip_fastq_in_close(quip_fastq_in_t*);
short_read_t* quip_fastq_read(quip_fastq_in_t*);


#endif
