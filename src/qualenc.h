/*
 * This file is part of quip.
 *
 * Copyright (c) 2011 by Daniel C. Jones <dcjones@cs.washington.edu>
 *
 */


/*
 * qual:
 * Statistical modeling of sequence quality scores.
 */

#ifndef QUIP_QUAL
#define QUIP_QUAL

#include "quip.h"
#include "parse.h"
#include <stdint.h>


typedef struct qualenc_t_ qualenc_t;

qualenc_t* qualenc_alloc(quip_block_writer_t writer, void* writer_data);
void       qualenc_free(qualenc_t*);

void qualenc_encode(qualenc_t*, const seq_t*);
void qualenc_clear(qualenc_t*);


#endif

