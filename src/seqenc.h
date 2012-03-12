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
#include "dist.h"
#include <stdlib.h>

typedef struct seqenc_t_ seqenc_t;

seqenc_t* seqenc_alloc_encoder(quip_writer_t writer, void* writer_data);
seqenc_t* seqenc_alloc_decoder(quip_reader_t writer, void* reader_data);
void      seqenc_free(seqenc_t*);

/* This is called to initialized the sequence motifs used when
 * encoding alignment. This must be called prior to any
 * calls to seqenc_encode_alignment. */
void seqenc_set_contigs(seqenc_t*, twobit_t** contigs, size_t n);

/* Update the contig sequences to the current maximum-likelihood
 * consensus sequence. */
void seqenc_get_contig_consensus(seqenc_t*, twobit_t** contigs);

void seqenc_encode_char_seq(seqenc_t*, const char*, size_t len);
void seqenc_encode_twobit_seq(seqenc_t*, const twobit_t*);

void seqenc_encode_alignment(
        seqenc_t* E,
        uint32_t contig_idx, uint32_t spos, uint8_t strand,
        const twobit_t* query);
void seqenc_flush(seqenc_t* E);

/* Optionally called to inform the decoder that the next n sequences should be
 * considered contigs. */
void seqenc_prepare_decoder(seqenc_t* E, uint32_t n, const uint32_t* lens);

void seqenc_decode(seqenc_t* E, seq_t* seq, size_t n);

void seqenc_start_decoder(seqenc_t* E);
void seqenc_reset_decoder(seqenc_t* E);


/* Number of mismatchisg reads before a contig position is flipped. */
extern const uint16_t mismatch_patch_factor;
extern const uint16_t mismatch_patch_cutoff;


#endif


