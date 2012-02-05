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

seqenc_t* seqenc_alloc_encoder(size_t k, quip_writer_t writer, void* writer_data);
seqenc_t* seqenc_alloc_decoder(size_t k, quip_reader_t writer, void* reader_data);
void      seqenc_free(seqenc_t*);

void seqenc_encode_char_seq(seqenc_t*, const char*, size_t len);
void seqenc_encode_twobit_seq(seqenc_t*, const twobit_t*);
void seqenc_encode_alignment(seqenc_t* E,
        size_t contig_idx, uint8_t strand,
        const sw_alignment_t* aln, const twobit_t* query);
void seqenc_flush(seqenc_t* E);

/* This is called optionally when we are expecting some number of 
 * of contig sequences. */
void seqenc_prepare_decoder(seqenc_t*, size_t contig_count);

bool seqenc_decode(seqenc_t* E, seq_t* seq, size_t n);
void seqenc_reset_decoder(seqenc_t* E);


#endif


