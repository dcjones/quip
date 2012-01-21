/*
 * This file is part of quip.
 *
 * Copyright (c) 2012 by Daniel C. Jones <dcjones@cs.washington.edu>
 *
 */


/*
 * ac:
 * A general purpose arithmetic encoder, following mostly from the
 * implementation described in "Introduction to Arithmetic Coding -- Theory and
 * Practice" by Amir Said.
 */


#ifndef QUIP_AC
#define QUIP_AC

#include "quip.h"
#include "dist.h"
#include <stdint.h>


typedef struct ac_t_ ac_t;

/* Allocate for encoding, decoding respectively. */
ac_t* ac_alloc_encoder(quip_writer_t writer, void* writer_data);
ac_t* ac_alloc_decoder(quip_reader_t reader, void* reader_data);
void  ac_free(ac_t*);

/* Encode symbol x from the given distribution */
void ac_encode(ac_t*, dist_t*, symb_t x);

/* Choose the final code value. */
void ac_flush_encoder(ac_t*);

/* Decode the next symbol with the given distribution. */
symb_t ac_decode(ac_t*, dist_t*);

#endif

