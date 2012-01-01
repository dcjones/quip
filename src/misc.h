/*
 * This file is part of quip.
 *
 * Copyright (c) 2011 by Daniel C. Jones <dcjones@cs.washington.edu>
 *
 * misc :
 * A few common functions, primarily for crashing whilst retaining our dignity.
 *
 */

#ifndef QUIP_MISC
#define QUIP_MISC

#include <stdio.h>

void or_die(int b, const char* msg);

void* malloc_or_die(size_t);
void* realloc_or_die(void*, size_t);
FILE* fopen_or_die(const char*, const char*);

#endif

