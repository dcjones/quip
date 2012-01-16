/*
 * This file is part of quip.
 *
 * Copyright (c) 2011 by Daniel C. Jones <dcjones@cs.washington.edu>
 *
 */


/*
 * assembler:
 * A dumb approximate de-novo assembler.
 *
 */


#ifndef QUIP_ASSEMBLER
#define QUIP_ASSEMBLER

#include "quip.h"
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>


typedef struct assembler_t_ assembler_t;

assembler_t* assembler_alloc(
        quip_block_writer_t writer, void* writer_data,
        size_t assemble_k, size_t align_k);
void         assembler_free(assembler_t*);

void assembler_add_seq(assembler_t*, const char* seq, size_t seqlen);
void assembler_assemble(assembler_t* A);

#endif

