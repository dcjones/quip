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
 * I do some truely unspeakable things with the proprocessor in order to
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
#define dec_size  12
#define dec_shift 12

#include "dist_template_on.h"
#include "dist_impl.c"
#include "dist_template_off.h"


#define DISTSIZE 4
#define dist_t dist4_t
#define cond_dist_t cond_dist4_t
#define dec_size  12
#define dec_shift 12

#include "dist_template_on.h"
#include "dist_impl.c"
#include "dist_template_off.h"


#define DISTSIZE 5
#define dist_t dist5_t
#define cond_dist_t cond_dist5_t
#define dec_size  12
#define dec_shift 12

#include "dist_template_on.h"
#include "dist_impl.c"
#include "dist_template_off.h"


#define DISTSIZE 72
#define dist_t dist72_t
#define cond_dist_t cond_dist72_t
#define dec_size  36
#define dec_shift 10

#include "dist_template_on.h"
#include "dist_impl.c"
#include "dist_template_off.h"


#define DISTSIZE 100
#define dist_t dist100_t
#define cond_dist_t cond_dist100_t
#define dec_size  36
#define dec_shift 10

#include "dist_template_on.h"
#include "dist_impl.c"
#include "dist_template_off.h"


