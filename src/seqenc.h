/*
 * This file is part of quip.
 *
 * Copyright (c) 2012 by Daniel C. Jones <dcjones@cs.washington.edu>
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

seqenc_t* seqenc_alloc_encoder(quip_writer_t writer, void* writer_data, const seqmap_t* ref);
seqenc_t* seqenc_alloc_decoder(quip_reader_t writer, void* reader_data, const seqmap_t* ref);
void      seqenc_free(seqenc_t*);

/* This is called to initialized the sequence motifs used when
 * encoding alignment. This must be called prior to any
 * calls to seqenc_encode_alignment. */
void seqenc_set_supercontig(seqenc_t*, const twobit_t* supercontig);

/* Update the contig sequences to the current maximum-likelihood
 * consensus sequence. */
void seqenc_get_supercontig_consensus(seqenc_t*, twobit_t* supercontig);

/* Encode/decode additional members of short_read. */
void seqenc_encode_extras(seqenc_t* E, const short_read_t* x, uint8_t quip_version);
void seqenc_decode_extras(seqenc_t* E, short_read_t* x, size_t seqlen, uint8_t quip_version);

void seqenc_encode_char_seq(seqenc_t*, const uint8_t*, size_t len);
void seqenc_encode_twobit_seq(seqenc_t*, const unsigned char* seq_str, const twobit_t* seq);

void seqenc_encode_reference_alignment(
        seqenc_t* E,
        const short_read_t*);

void seqenc_encode_alignment(
        seqenc_t* E,
        uint32_t spos, uint8_t strand,
        const unsigned char *query_str,
        const twobit_t* query);

size_t seqenc_finish(seqenc_t* E);
void   seqenc_flush(seqenc_t* E);

void seqenc_decode(seqenc_t* E, short_read_t* seq, size_t n);

void seqenc_start_decoder(seqenc_t* E);
void seqenc_reset_decoder(seqenc_t* E);


/* Number of mismatchisg reads before a contig position is flipped. */
extern const uint16_t mismatch_patch_factor;
extern const uint16_t mismatch_patch_cutoff;


#endif


