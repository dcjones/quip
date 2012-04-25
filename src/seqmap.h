/*
 * This file is part of quip.
 *
 * Copyright (c) 2011 by Daniel C. Jones <dcjones@cs.washington.edu>
 *
 */

/*
 * seqmap:
 * Map sequence names to sequences, maintaining a set of named sequences used
 * in reference based compression of aligned reads.
 *
 */

#ifndef QUIP_SEQMAP
#define QUIP_SEQMAP

#include "twobit.h"
#include <stdio.h>

typedef struct seqmap_t_ seqmap_t;

seqmap_t* seqmap_alloc();
void      seqmap_clear(seqmap_t*);
void      seqmap_free(seqmap_t*);
void      seqmap_read_fasta(seqmap_t*, FILE*);
const twobit_t* seqmap_get(seqmap_t*, const char* seqname);

#endif
