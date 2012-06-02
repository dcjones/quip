/*
 * This file is part of quip.
 *
 * Copyright (c) 2012 by Daniel C. Jones <dcjones@cs.washington.edu>
 *
 */


/*
 * samoptenc:
 * Compression and decompression of optional SAM fields.
 *
 * Note: definitions are in samopt.c
 */


#ifndef QUIP_SAMOPTENC
#define QUIP_SAMOPTENC

#include "quip.h"

typedef struct samoptenc_t_ samoptenc_t;

samoptenc_t* samoptenc_alloc_encoder(quip_writer_t writer, void* writer_data);
samoptenc_t* samoptenc_alloc_decoder(quip_reader_t reader, void* reader_data);
void samoptenc_free(samoptenc_t* E);

void samoptenc_encode(samoptenc_t* E, const samopt_table_t* T);
void samoptenc_decode(samoptenc_t* E, samopt_table_t* T);

size_t samoptenc_finish(samoptenc_t* E);
void   samoptenc_flush(samoptenc_t* E);

void samoptenc_start_decoder(samoptenc_t* E);
void samoptenc_reset_decoder(samoptenc_t* E);

/* Update the checksum with the last encoded aux data. */
uint64_t samoptenc_crc64_update(const samoptenc_t* E, uint64_t crc);

#endif
