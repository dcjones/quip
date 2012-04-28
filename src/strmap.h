/*
 * This file is part of quip.
 *
 * Copyright (c) 2012 by Daniel C. Jones <dcjones@cs.washington.edu>
 *
 */


/*
 * strmap:
 * map string to sequential integers
 */

#ifndef QUIP_STRMAP
#define QUIP_STRMAP

#include "quip.h"

typedef struct strmap_t_ strmap_t;

strmap_t* strmap_alloc();
void      strmap_free(strmap_t*);
size_t    strmap_size(const strmap_t*);

/* get the integer mapped to the string, inserting it if it does not exist */
uint32_t  strmap_get(strmap_t*, const str_t*);

#endif
