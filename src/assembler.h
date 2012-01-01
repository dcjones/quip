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

#include <stdlib.h>
#include <stdio.h>


typedef struct assembler_t_ assembler_t;

assembler_t* assembler_alloc(size_t assemble_k, size_t align_k);
void         assembler_free(assembler_t*);

void assembler_add_seq(assembler_t*, const char* seq, size_t seqlen);
void assembler_write(assembler_t*, FILE*);

#endif

