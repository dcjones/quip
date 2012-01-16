/*
 * This file is part of quip.
 *
 * Copyright (c) 2011 by Daniel C. Jones <dcjones@cs.washington.edu>
 *
 */


/*
 * seqenc:
 * Compression of nucleotide sequences.
 */

#ifndef QUIP_SEQENC
#define QUIP_SEQENC

#include "quip.h"
#include "twobit.h"
#include "sw.h"
#include <stdlib.h>

typedef struct seqenc_t_ seqenc_t;

seqenc_t* seqenc_alloc(size_t k, quip_block_writer_t writer, void* writer_data);
void      seqenc_free(seqenc_t*);

void seqenc_encode_char_seq(seqenc_t*, const char*);
void seqenc_encode_twobit_seq(seqenc_t*, const twobit_t*);
void seqenc_encode_alignment(seqenc_t* E,
        size_t contig_idx, uint8_t strand,
        const sw_alignment_t* aln, const twobit_t* query);

void seqenc_flush(seqenc_t* E);


#endif

