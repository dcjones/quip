/*
 * This file is part of quip.
 *
 * Copyright (c) 2012 by Daniel C. Jones <dcjones@cs.washington.edu>
 *
 */


/*
 * ac:
 * A general purpose arithmetic encoder, following mostly from the
 * implementation described in "Introduction to Arithmetic Coding -- Theory and
 * Practice" by Amir Said.
 */


#ifndef QUIP_AC
#define QUIP_AC

#include "quip.h"
#include <stdint.h>

typedef struct ac_t_
{
    /* coder state */
    uint32_t b; /* base */
    uint32_t l; /* length */
    uint32_t v; /* value */

    /* input or output buffer (depending on whether we are decoding or encoding,
     * respectively). */
    uint8_t* buf;

    /* size allocated to buf */
    size_t buflen;

    /* index next vacant buffer position */
    size_t bufpos;

    /* available indput (for decoding) */
    size_t bufavail;

    /* callback function for encoder output */
    quip_writer_t writer;
    void* writer_data;

    /* callback function for decoder input */
    quip_reader_t reader;
    void* reader_data;
} ac_t;

extern const uint32_t min_length;
extern const uint32_t max_length;

/* Allocate for encoding, decoding respectively. */
ac_t* ac_alloc_encoder(quip_writer_t writer, void* writer_data);
ac_t* ac_alloc_decoder(quip_reader_t reader, void* reader_data);
void  ac_free(ac_t*);

/* Choose the final code value and return the number of compressed bytes. */
size_t ac_finish_encoder(ac_t*);

/* Write buffered compressed bases and reset the encoder. */
void ac_flush_encoder(ac_t*);

/* This must be called following a call to alloc_decoder or
 * reset_decoder prior to  any symbols are decoded. */
void ac_start_decoder(ac_t*);

/* Start over decoding. */
void ac_reset_decoder(ac_t*);

/* If you don't know what these do, don't call them. */
void ac_renormalize_encoder(ac_t*);
void ac_renormalize_decoder(ac_t*);
void ac_propogate_carry(ac_t*);

#endif

