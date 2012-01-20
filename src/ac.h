/*
 * This file is part of quip.
 *
 * Copyright (c) 2011 by Daniel C. Jones <dcjones@cs.washington.edu>
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
#include "cumdist.h"
#include <stdint.h>


/* Encoder */

typedef struct ac_t_ ac_t;

ac_t* ac_alloc(quip_block_writer_t writer, void* writer_data);
void  ac_free(ac_t*);

/* Update the arithmetic coder with the symbol having probability p and
 * cumulative probability P. */
void ac_update(ac_t*, uint32_t p, uint32_t P);

/* Choose the final code value. */
void ac_flush(ac_t*);



/* Decoder */

typedef struct dec_t_ dec_t;

dec_t* dec_alloc(quip_reader_t reader, void* reader_data);
void dec_free(dec_t*);

size_t dec_next(dec_t*, cumdist_t* C);

#endif

