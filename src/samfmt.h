/*
 * This file is part of quip.
 *
 * Copyright (c) 2012 by Daniel C. Jones <dcjones@cs.washington.edu>
 *
 */

#ifndef QUIP_SAMFMT
#define QUIP_SAMFMT

#include "quip.h"

typedef struct quip_sam_out_t_ quip_sam_out_t;

quip_sam_out_t* quip_sam_out_open(
                    quip_writer_t writer,
                    void*         writer_data,
                    quip_opt_t    opts);

void quip_sam_out_close(quip_sam_out_t*);
void quip_sam_write(quip_sam_out_t*, short_read_t*);


typedef struct quip_sam_in_t_ quip_sam_in_t;

quip_sam_in_t* quip_sam_in_open(
                    quip_reader_t reader,
                    void*         reader_data,
                    quip_opt_t opts);

void quip_sam_in_close(quip_sam_in_t*);
short_read_t* quip_sam_read(quip_sam_in_t*);


#endif

