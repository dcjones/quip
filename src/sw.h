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


#ifndef QUIP_SW_H
#define QUIP_SW_H

#include "twobit.h"

typedef struct sw_t_ sw_t;


typedef enum edit_op_t_
{
    EDIT_MATCH = 0,
    EDIT_Q_GAP,
    EDIT_S_GAP,
    EDIT_MISMATCH
} edit_op_t;

typedef struct sw_alignment_t_
{
    size_t len;  /* length of the stored alignment */
    size_t size; /* space allocated in ops */
    int spos;    /* position within the subject that the alignment begins */
    edit_op_t* ops;
} sw_alignment_t;

void print_alignment(FILE*, const sw_alignment_t*);

sw_t* sw_alloc(const twobit_t* subject);
void  sw_set_subject(sw_t*, const twobit_t* subject);
void  sw_free(sw_t*);

int sw_seeded_align(sw_t* sw, const twobit_t* query,
                    int spos, int qpos, int seedlen);


/* Store the alignment found in the last call to sw_seeded_align. */
void sw_trace(sw_t* sw, sw_alignment_t* aln);

/* Band size used by the striped smith-waterman approach */
extern const int sw_band_width;

#endif

