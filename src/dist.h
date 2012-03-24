

#ifndef QUIP_DIST
#define QUIP_DIST

#include <stdlib.h>
#include <stdint.h>

extern const size_t dist_length_shift;

/* Integer used to represent symbols in the alphabet. */
typedef uint32_t symb_t;

#define DISTSIZE 2
#define dist_t dist2_t
#define cond_dist_t cond_dist2_t

#include "dist_template_on.h"
#include "dist_impl.h"
#include "dist_template_off.h"


#define DISTSIZE 4
#define dist_t dist4_t
#define cond_dist_t cond_dist4_t

#include "dist_template_on.h"
#include "dist_impl.h"
#include "dist_template_off.h"


#define DISTSIZE 16 
#define dist_t dist16_t
#define cond_dist_t cond_dist16_t

#include "dist_template_on.h"
#include "dist_impl.h"
#include "dist_template_off.h"


#define DISTSIZE 41 
#define dist_t dist41_t
#define cond_dist_t cond_dist41_t

#include "dist_template_on.h"
#include "dist_impl.h"
#include "dist_template_off.h"


#define DISTSIZE 50
#define dist_t dist50_t
#define cond_dist_t cond_dist50_t

#include "dist_template_on.h"
#include "dist_impl.h"
#include "dist_template_off.h"


#define DISTSIZE 64
#define dist_t dist64_t
#define cond_dist_t cond_dist64_t

#include "dist_template_on.h"
#include "dist_impl.h"
#include "dist_template_off.h"


#define DISTSIZE 128
#define dist_t dist128_t
#define cond_dist_t cond_dist128_t

#include "dist_template_on.h"
#include "dist_impl.h"
#include "dist_template_off.h"


#define DISTSIZE 256
#define dist_t dist256_t
#define cond_dist_t cond_dist256_t

#include "dist_template_on.h"
#include "dist_impl.h"
#include "dist_template_off.h"

void dist_encode_uint32(ac_t* ac, cond_dist256_t* d, uint32_t x);
uint32_t dist_decode_uint32(ac_t* ac, cond_dist256_t* d);

void dist_encode_uint64(ac_t* ac, cond_dist256_t* d, uint64_t x);
uint64_t dist_decode_uint64(ac_t* ac, cond_dist256_t* d);

#endif

