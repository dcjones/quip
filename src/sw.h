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
    EDIT_MATCH,
    EDIT_MISMATCH,
    EDIT_Q_GAP,
    EDIT_S_GAP
} edit_op_t;

typedef struct sw_alignment_t_
{
    size_t len;  /* length of the stored alignment */
    size_t size; /* space allocated in ops */
    int spos; /* position within the subject that the alignment begins */
    edit_op_t* ops;
} sw_alignment_t;

sw_t* sw_alloc(const twobit_t* subject);
void  sw_free(sw_t*);

int sw_seeded_align(sw_t* sw, const twobit_t* query,
                    int spos, int qpos, int seedlen);


/* Store the alignment found in the last call to sw_seeded_align. */
void sw_trace(sw_t* sw, sw_alignment_t* aln);


/* align a query sequence to the subject, with the given anchor */
//int sw_align(sw_t* subject, const twobit_t* query,
             //int spos, int qpos, int seedlen);


/* TODO: a function to align not-twobit sequences */

/* TODO: a function to actually retrieve the alignment */



/* these functions are for testing */
//sw_t* sw_alloc_char(const char* subject);
//int sw_align_char(sw_t* subject, const char* query);

#endif

