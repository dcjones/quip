/*
 * This file is part of quip.
 *
 * Copyright (c) 2012 by Daniel C. Jones <dcjones@cs.washington.edu>
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
        quip_writer_t   writer,
        void*           writer_data,
        bool            quick,
        const seqmap_t* ref);

void assembler_free(assembler_t*);

void assembler_clear_contigs(assembler_t*);

void   assembler_add_seq(assembler_t*, const short_read_t* seq);
size_t assembler_finish(assembler_t* A);
void   assembler_flush(assembler_t* A);


/* disassemble */
typedef struct disassembler_t_ disassembler_t;

disassembler_t* disassembler_alloc(
    quip_reader_t reader,
    void* reader_data,
    bool quick,
    const seqmap_t* ref);

void disassembler_free(disassembler_t*);

void disassembler_read(disassembler_t*, short_read_t* x, size_t n);
void disassembler_reset(disassembler_t*);


#endif

