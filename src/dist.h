

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


#define DISTSIZE 8
#define dist_t dist8_t
#define cond_dist_t cond_dist8_t

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


#define DISTSIZE 100
#define dist_t dist100_t
#define cond_dist_t cond_dist100_t

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

#endif

