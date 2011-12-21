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

typedef struct qualmodel_t_ qualmodel_t;


qualmodel_t* qualmodel_alloc();
void         qualmodel_free(qualmodel_t*);

/* Update the quality model with the given n sequences. */
void qualmodel_update(qualmodel_t*, const seq_t* x);


#endif

