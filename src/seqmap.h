/*
 * This file is part of quip.
 *
 * Copyright (c) 2012 by Daniel C. Jones <dcjones@cs.washington.edu>
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

#include "quip.h"
#include "twobit.h"
#include <stdio.h>

seqmap_t* seqmap_alloc();
void      seqmap_clear(seqmap_t*);
void      seqmap_free(seqmap_t*);
void      seqmap_read_fasta(seqmap_t*, FILE*);
size_t    seqmap_size(const seqmap_t*);
const twobit_t* seqmap_get(const seqmap_t*, const char* seqname, size_t* idx);
uint64_t        seqmap_crc64(const seqmap_t*);

#endif
