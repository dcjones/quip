/*
 * This file is part of quip.
 *
 * Copyright (c) 2012 by Daniel C. Jones <dcjones@cs.washington.edu>
 *
 */


/*
 * dist:
 * Efficient representation of discrete probability distributions for arithmetic
 * coding.
 *
 * I do some truly unspeakable things with the preprocessor in order to
 * slightly improve the caching behavior of the arithmetic coder. This is
 * basically C++ template metagrogramming but in C. I use defines to declare
 * several distribution objects of fixed size.
 */



#include "dist.h"
#include "misc.h"
#include <string.h>
#include <assert.h>


/* The distribution is updated every "initial_update_factor * n" new
 * observations. */
static const size_t update_delay_factor = 1;
const size_t dist_length_shift = 15;
static const size_t max_count = 1 << 15;


#define DISTSIZE 2
#define dist_t dist2_t
#define cond_dist_t cond_dist2_t

#include "dist_template_on.h"
#include "dist_impl.c"
#include "dist_template_off.h"


#define DISTSIZE 4
#define dist_t dist4_t
#define cond_dist_t cond_dist4_t

#include "dist_template_on.h"
#include "dist_impl.c"
#include "dist_template_off.h"


#define DISTSIZE 16
#define dist_t dist16_t
#define cond_dist_t cond_dist16_t

#include "dist_template_on.h"
#include "dist_impl.c"
#include "dist_template_off.h"


#define DISTSIZE 50
#define dist_t dist50_t
#define cond_dist_t cond_dist50_t

#include "dist_template_on.h"
#include "dist_impl.c"
#include "dist_template_off.h"


#define DISTSIZE 64
#define dist_t dist64_t
#define cond_dist_t cond_dist64_t

#include "dist_template_on.h"
#include "dist_impl.c"
#include "dist_template_off.h"


#define DISTSIZE 128
#define dist_t dist128_t
#define cond_dist_t cond_dist128_t

#include "dist_template_on.h"
#include "dist_impl.c"
#include "dist_template_off.h"


#define DISTSIZE 256
#define dist_t dist256_t
#define cond_dist_t cond_dist256_t

#include "dist_template_on.h"
#include "dist_impl.c"
#include "dist_template_off.h"


void uint32_enc_init(uint32_enc_t* d)
{
    cond_dist256_init(d, 5 * 256);
}

void uint32_enc_free(uint32_enc_t* d)
{
    cond_dist256_free(d);
}

void uint32_enc_encode(ac_t* ac, uint32_enc_t* d, uint32_t x)
{
    uint32_t y, y0;

    y0 = 0;
    y = x >> 24;
    cond_dist256_encode(ac, d, y0, y);

    y0 = y;
    y = (x >> 16) & 0xff;
    cond_dist256_encode(ac, d, 0x100 | y0, y);

    y0 = y;
    y = (x >> 8) & 0xff;
    cond_dist256_encode(ac, d, 0x200 | y0, y);

    y0 = y;
    y = x & 0xff;
    cond_dist256_encode(ac, d, 0x300 | y0, y);
}

uint32_t uint32_enc_decode(ac_t* ac, uint32_enc_t* d)
{
    uint32_t y, y0, x;

    y0 = 0;
    x = cond_dist256_decode(ac, d, y0);

    y0 = x;
    y = cond_dist256_decode(ac, d, 0x100 | y0);
    x = (x << 8) | y;

    y0 = y;
    y = cond_dist256_decode(ac, d, 0x200 | y0);
    x = (x << 8) | y;

    y0 = y;
    y = cond_dist256_decode(ac, d, 0x300 | y0);
    x = (x << 8) | y;

    return x;
}

void uint64_enc_init(uint64_enc_t* d)
{
    cond_dist256_init(d, 9 * 256);
}

void uint64_enc_free(uint64_enc_t* d)
{
    cond_dist256_free(d);
}

void uint64_enc_encode(ac_t* ac, uint64_enc_t* d, uint64_t x)
{
    uint64_t y0 = 0, y;
    size_t k = 0;
    while (x >= 0x80) {
        y  = (uint8_t) x | 0x80;
        cond_dist256_encode(ac, d, (k << 8) | y0, y);
        x >>= 7;
        ++k;
        y0 = y;
    }

    cond_dist256_encode(ac, d, (k << 8) | y0, x);
}


uint64_t uint64_enc_decode(ac_t* ac, uint64_enc_t* d)
{
    uint64_t x, y0;

    x = 0;
    y0 = 0;
    size_t k = 0;
    do {
        y0 = cond_dist256_decode(ac, d, (k << 8) | y0);
        x |= (y0 & 0x7f) << (k * 7);
        ++k;
    } while (y0 & 0x80);

    return x;
}

