

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

typedef cond_dist256_t uint32_enc_t;
void uint32_enc_init(uint32_enc_t*);
void uint32_enc_free(uint32_enc_t*);
void uint32_enc_encode(ac_t*, uint32_enc_t*, uint32_t);
uint32_t uint32_enc_decode(ac_t*, uint32_enc_t*);

typedef cond_dist256_t uint64_enc_t;
void uint64_enc_init(uint64_enc_t*);
void uint64_enc_free(uint64_enc_t*);
void uint64_enc_encode(ac_t*, uint64_enc_t*, uint64_t);
uint64_t uint64_enc_decode(ac_t*, uint64_enc_t*);

#endif

