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

#include "config.h"
#include <stdio.h>
#include <stdint.h>

void or_die(int b, const char* msg);

void* malloc_or_die(size_t);
void* realloc_or_die(void*, size_t);
FILE* fopen_or_die(const char*, const char*);

uint32_t strhash(const char*, size_t len);

#if HAVE_PREFETCH
#define prefetch(p) __builtin_prefetch(p)
#else
#define prefetch(p)
#endif

#endif

