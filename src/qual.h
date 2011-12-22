/*
 * This file is part of fazer.
 *
 * Copyright (c) 2011 by Daniel C. Jones <dcjones@cs.washington.edu>
 *
 */


/*
 * qual:
 * Statistical modeling of sequence quality scores.
 */

#ifndef FAZER_QUAL
#define FAZER_QUAL

#include "parse.h"
#include <stdint.h>


typedef struct qualenc_t_ qualenc_t;

typedef void (*qualenc_writer_t) (void*, uint8_t*, size_t);

qualenc_t* qualenc_alloc(qualenc_writer_t, void* writer_data);
void       qualenc_free(qualenc_t*);

void qualenc_encode(qualenc_t*, const seq_t*);
void qualenc_finish(qualenc_t*);


#endif

