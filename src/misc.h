/*
 * This file is part of quip.
 *
 * Copyright (c) 2012 by Daniel C. Jones <dcjones@cs.washington.edu>
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
#include <stdarg.h>

void* malloc_or_die(size_t);
void* realloc_or_die(void*, size_t);
FILE* fopen_or_die(const char*, const char*);

#if HAVE_PREFETCH
#define prefetch(p, rw, locality) __builtin_prefetch(p, rw, locality)
#else
#define prefetch(p, rw, locality)
#endif

#ifndef __attribute__
#if __GNUC__ < 2 || (__GNUC__ == 2 && __GNUC_MINOR__ < 8) || __STRICT_ANSI__
#define __attribute__(x)
#endif
#endif

#define ATTRIB_UNUSED __attribute__((__unused__))


/* Windows reads/writes in "text mode" by default. This is confusing 
 * and wrong, so we need to disable it.
 */
#if defined(MSDOS) || defined(OS2) || defined(WIN32) || defined(__CYGWIN__)
#  include <fcntl.h>
#  include <io.h>
#  define SET_BINARY_MODE(file) setmode(fileno(file), O_BINARY)
#else
#  define SET_BINARY_MODE(file)
#endif


/* Generic hashing, using MurmurHash3 */
uint32_t murmurhash3(const uint8_t* data, size_t len);

unsigned char complement(unsigned char c);
void str_revcomp(unsigned char* seq, size_t n);
void str_rev(unsigned char* seq, size_t n);

#ifndef HAVE_VASPRINTF
int vasprintf(char **ret, const char *format, va_list args);
#endif

#ifndef HAVE_VASPRINTF
int asprintf(char **ret, const char *format, ...);
#endif

#endif

