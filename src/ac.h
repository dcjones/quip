/*
 * This file is part of quip.
 *
 * Copyright (c) 2011 by Daniel C. Jones <dcjones@cs.washington.edu>
 *
 */


/*
 * ac:
 * A general purpose arithmetic encoder.
 */


#ifndef QUIP_AC
#define QUIP_AC

#include "quip.h"
#include <stdint.h>


typedef struct ac_t_ ac_t;

ac_t* ac_alloc(quip_block_writer_t writer, void* writer_data);
void  ac_free(ac_t*);

/* Update the arithmetic coder with the symbol having probability p and
 * cumulative probability P. */
void ac_update(ac_t*, uint32_t p, uint32_t P);

#endif

