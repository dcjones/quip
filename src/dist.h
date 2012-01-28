

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
#define dec_size  12
#define dec_shift 12

#include "dist_template_on.h"
#include "dist_impl.h"
#include "dist_template_off.h"


#define DISTSIZE 4
#define dist_t dist4_t
#define cond_dist_t cond_dist4_t
#define dec_size  12
#define dec_shift 12

#include "dist_template_on.h"
#include "dist_impl.h"
#include "dist_template_off.h"


#define DISTSIZE 8
#define dist_t dist8_t
#define cond_dist_t cond_dist8_t
#define dec_size  12
#define dec_shift 12

#include "dist_template_on.h"
#include "dist_impl.h"
#include "dist_template_off.h"


#define DISTSIZE 72
#define dist_t dist72_t
#define cond_dist_t cond_dist72_t
#define dec_size  36
#define dec_shift 10

#include "dist_template_on.h"
#include "dist_impl.h"
#include "dist_template_off.h"


#define DISTSIZE 100
#define dist_t dist100_t
#define cond_dist_t cond_dist100_t
#define dec_size  36
#define dec_shift 10

#include "dist_template_on.h"
#include "dist_impl.h"
#include "dist_template_off.h"


#define DISTSIZE 125
#define dist_t dist125_t
#define cond_dist_t cond_dist125_t
#define dec_size  12
#define dec_shift 12

#include "dist_template_on.h"
#include "dist_impl.h"
#include "dist_template_off.h"


#define DISTSIZE 128
#define dist_t dist128_t
#define cond_dist_t cond_dist128_t
#define dec_size  36
#define dec_shift 10

#include "dist_template_on.h"
#include "dist_impl.h"
#include "dist_template_off.h"


#define DISTSIZE 256
#define dist_t dist256_t
#define cond_dist_t cond_dist256_t
#define dec_size  68
#define dec_shift 9

#include "dist_template_on.h"
#include "dist_impl.h"
#include "dist_template_off.h"

#endif

