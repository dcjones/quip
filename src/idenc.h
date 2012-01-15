/*
 * This file is part of quip.
 *
 * Copyright (c) 2011 by Daniel C. Jones <dcjones@cs.washington.edu>
 *
 */


/*
 * idenc:
 * Compression and decompression of sequence ids.
 */

#ifndef QUIP_IDENC
#define QUIP_IDENC

#include "quip.h"
#include "parse.h"
#include <stdint.h>


typedef struct idenc_t_ idenc_t;

idenc_t* idenc_alloc(quip_block_writer_t writer, void* writer_data);
void     idenc_free(idenc_t*);

void idenc_encode(idenc_t*, const seq_t*);
void idenc_flush(idenc_t*);

#endif

