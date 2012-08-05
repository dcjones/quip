/*
 * This file is part of quip.
 *
 * Copyright (c) 2012 by Daniel C. Jones <dcjones@cs.washington.edu>
 *
 */


/*
 * samopt_table functions that are not exposed as part of the external API.
 */

#ifndef QUIP_SAMOPT
#define QUIP_SAMOPT

#include "quip.h"
#include "sam/bam.h"

void samopt_table_bam_dump(const samopt_table_t*, bam1_t*);

#endif

