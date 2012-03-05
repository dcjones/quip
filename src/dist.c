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

/* Code for computing dec_size and dec_shift:
 *
 *      size_t dec_bits = 3;
 *      while (DISTSIZE > (1U << (dec_bits + 2))) ++dec_bits;
 *      dec_size  = (1 << dec_bits) + 4;
 *      dec_shift = dist_length_shift - dec_bits;
 */


#define DISTSIZE 2
#define dist_t dist2_t
#define cond_dist_t cond_dist2_t

#include "dist_template_on.h"
#include "dist_impl.c"
#include "dist_template_off.h"

 
#define DISTSIZE 3
#define dist_t dist3_t
#define cond_dist_t cond_dist3_t

#include "dist_template_on.h"
#include "dist_impl.c"
#include "dist_template_off.h"


#define DISTSIZE 4
#define dist_t dist4_t
#define cond_dist_t cond_dist4_t

#include "dist_template_on.h"
#include "dist_impl.c"
#include "dist_template_off.h"


#define DISTSIZE 8
#define dist_t dist8_t
#define cond_dist_t cond_dist8_t

#include "dist_template_on.h"
#include "dist_impl.c"
#include "dist_template_off.h"


#define DISTSIZE 16
#define dist_t dist16_t
#define cond_dist_t cond_dist16_t

#include "dist_template_on.h"
#include "dist_impl.c"
#include "dist_template_off.h"


#define DISTSIZE 41
#define dist_t dist41_t
#define cond_dist_t cond_dist41_t

#include "dist_template_on.h"
#include "dist_impl.c"
#include "dist_template_off.h"


#define DISTSIZE 50
#define dist_t dist50_t
#define cond_dist_t cond_dist50_t

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


/* Functions for general-purpose markov-chain
 * encoding of unsigned integers.
 */

void dist_encode_uint32(ac_t* ac, cond_dist256_t* d, uint32_t x)
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

uint32_t dist_decode_uint32(ac_t* ac, cond_dist256_t* d)
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


void dist_encode_uint64(ac_t* ac, cond_dist256_t* d, uint64_t x)
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


uint64_t dist_decode_uint64(ac_t* ac, cond_dist256_t* d)
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

