/*
 * This file is part of quip.
 *
 * Copyright (c) 2012 by Daniel C. Jones <dcjones@cs.washington.edu>
 *
 * parse :
 * A parser for FASTQ files.
 *
 */

#ifndef QUIP_PARSE_H
#define QUIP_PARSE_H

#include <stdio.h>


typedef struct
{
    char*  s;    /* null-terminated string */
    size_t n;    /* length of s */
    size_t size; /* bytes allocated for s */
} str_t;

void fastq_expand_str(str_t* s);
void fastq_reserve_str(str_t* s, size_t n);


typedef struct
{
    str_t id1;
    str_t seq;
    str_t id2;
    str_t qual;
} seq_t;


seq_t* fastq_alloc_seq();
void fastq_copy_seq(seq_t* dest, const seq_t* src);
void fastq_free_seq(seq_t*);


typedef struct
{
    FILE*  file;
    int    state;
    char*  buf;
    char*  c;
} fastq_t;


fastq_t* fastq_open(FILE*);
void fastq_close(fastq_t*);
int  fastq_next(fastq_t*, seq_t*);
void fastq_rewind(fastq_t*);

#endif

