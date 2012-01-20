/*
 * This file is part of quip.
 *
 * Copyright (c) 2012 by Daniel C. Jones <dcjones@cs.washington.edu>
 *
 */

#ifndef QUIP_QUIP
#define QUIP_QUIP

#include "parse.h"
#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>

typedef void   (*quip_block_writer_t) (void*, const uint8_t*, size_t);
typedef size_t (*quip_reader_t) (void*, uint8_t*, size_t);

typedef struct quip_compressor_t_ quip_compressor_t;
quip_compressor_t* quip_comp_alloc(quip_block_writer_t, void* writer_data, bool quick);
void               quip_comp_addseq(quip_compressor_t*, seq_t*);
void               quip_comp_flush(quip_compressor_t*);
void               quip_comp_free(quip_compressor_t*);

void quip_write_header(FILE*);

extern bool verbose;

#endif


