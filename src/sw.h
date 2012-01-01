/*
 * This file is part of quip.
 *
 * Copyright (c) 2011 by Daniel C. Jones <dcjones@cs.washington.edu>
 *
 */

/*
 * sw :
 * Local alignments of nucleotide sequences via Smith-Waterman.
 *
 */


/* Here's what I want to be able to do:
 *  Allocate a sw_t for each contig.
 *
 *  For this to work, we need to add contstraints when performing the alignment.
 *  In particular, we need to be able to specify a seed, and bound the subject
 *  (to disallow very long gaps).
 *
 */


#ifndef QUIP_SW_H
#define QUIP_SW_H

#include "twobit.h"

typedef struct sw_t_ sw_t;

sw_t* sw_alloc(const twobit_t* subject);
void  sw_free(sw_t*);

/* align a query sequence to the subject, with the given anchor */
int sw_align(sw_t* subject, const twobit_t* query,
             int spos, int qpos, int seedlen);


/* TODO: a function to align not-twobit sequences */

/* TODO: a function to actually retrieve the alignment */



/* these functions are for testing */
sw_t* sw_alloc_char(const char* subject);
int sw_align_char(sw_t* subject, const char* query);

#endif

