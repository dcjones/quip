/*
 * This file is part of quip.
 *
 * Copyright (c) 2011 by Daniel C. Jones <dcjones@cs.washington.edu>
 *
 */


#include "sw.h"
#include "misc.h"
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <limits.h>


void print_alignment(FILE* fout, const sw_alignment_t* aln)
{
    size_t i;
    for (i = 0; i < aln->len; ++i) {
        switch (aln->ops[i]) {
            case EDIT_MATCH:    fputc('M', fout); break;
            case EDIT_MISMATCH: fputc('N', fout); break;
            case EDIT_Q_GAP:    fputc('Q', fout); break;
            case EDIT_S_GAP:    fputc('S', fout); break;
        }
    }
    fputc('\n', fout);
}



typedef uint16_t score_t;


struct sw_t_
{
    uint8_t* subject;
    uint8_t* query;

    int n; /* subject length */
    int m; /* query length */

    int qlen;
    int spos;

    int last_spos;

    /* start within the subject sequence of the alignment */

    /* overall minimum cost alignment */
    score_t* F;
};


/* Scoring scheme */
static const score_t score_inf      = UINT16_MAX / 2;
static const score_t score_q_gap    = 2;
static const score_t score_s_gap    = 4;
static const score_t score_match    = 1;
static const score_t score_mismatch = 3;

static const int band_width = 1;
static const int colsize = 3; /* 1 + 2 * band_width */



sw_t* sw_alloc(const twobit_t* subject)
{
    sw_t* sw = malloc_or_die(sizeof(sw_t));

    sw->n = (int) twobit_len(subject);

    sw->subject = malloc_or_die((sw->n + 2 * band_width) * sizeof(uint8_t));

    int i;
    /* the lead of the subject sequence is padded to slightly simplify the
     * alignment */
    for (i = 0; i < band_width; ++i) sw->subject[i] = 5; /* '5' to ensure mismatches */
    for (i = 0; i < sw->n; ++i) sw->subject[band_width + i] = twobit_get(subject, i);
    for (i = 0; i < band_width; ++i) sw->subject[sw->n + i] = 5;

    /* We don't know the query length in advance, so these are all allocated
     * and/or expanded when needed. */
    sw->m = 0;
    sw->query = NULL;
    sw->F = NULL;

    return sw;
}


void sw_set_subject(sw_t* sw, const twobit_t* subject)
{
    int len = twobit_len(subject);

    if (len != sw->n) {
        free(sw->subject);
        sw->subject = malloc_or_die((len + band_width) * sizeof(uint8_t));
        sw->n = len;
    }

    int i;
    /* the lead of the subject sequence is padded to slightly simplify the
     * alignment */
    for (i = 0; i < band_width; ++i) sw->subject[i] = 5; /* '5' to ensure mismatches */
    for (i = 0; i < sw->n; ++i) sw->subject[band_width + i] = twobit_get(subject, i);
}


void sw_free(sw_t* sw)
{
    free(sw->subject);
    free(sw->query);
    free(sw->F);
    free(sw);
}


static inline score_t min_score(score_t a, score_t b)
{
    return a < b ? a : b;
}


/* Align the query block between qu and qv, inclusively. */ 
static void sw_align_sub(sw_t* sw, int qu, int qv)
{
    int colstart, colend;
    int i, j, k;

    for (i = qu; i <= qv; ++i) {
        j = colstart = (i + 1) * colsize;
        colend = colstart + colsize;

        /* special case for the first entry in the column */

        sw->F[j] = sw->F[j - colsize] +
            (sw->subject[sw->spos + i] == sw->query[i] ?
             score_match : score_mismatch);

        sw->F[j] = min_score(sw->F[j], sw->F[j - colsize + 1] + score_s_gap);

        ++j;


        /* middle cells */

        for (k = 1; j < colend - 1; ++j, ++k) {
            sw->F[j] = sw->F[j - colsize] +
                (sw->subject[sw->spos + i + k] == sw->query[i] ?
                 score_match : score_mismatch);

            sw->F[j] = min_score(sw->F[j], sw->F[j - colsize + 1] + score_s_gap);

            sw->F[j] = min_score(sw->F[j], sw->F[j - 1] + score_q_gap);
        }


        /* special case for the last entry in the column */

        sw->F[j] = sw->F[j - colsize] +
            (sw->subject[sw->spos + i + k] == sw->query[i] ?
             score_match : score_mismatch);

        sw->F[j] = min_score(sw->F[j], sw->F[j - 1] + score_q_gap);
    }
}


static void sw_ensure_query_len(sw_t* sw, int m)
{
    if (sw->m < m)  {
        sw->m = m;
        sw->query = realloc_or_die(sw->query, sw->m * sizeof(uint8_t));

        free(sw->F);
        sw->F = malloc_or_die((m + 1) * colsize * sizeof(score_t));
    }
}




int sw_seeded_align(sw_t* sw, const twobit_t* query,
                    int spos, int qpos, int seedlen)
{
    int qlen = sw->qlen = (int) twobit_len(query);

    if (qpos > spos || (qlen - qpos) > (sw->n - spos)) return score_inf;


    sw_ensure_query_len(sw, qlen);
    sw->spos = spos - qpos;
    int i, j;

    for (i = 0; i < qlen; ++i) sw->query[i] = twobit_get(query, i);

    /* align up to the seed */
    memset(sw->F, 0, colsize * sizeof(score_t));

    /* account for the case in which we are aligning right at the start of the
     * contig */
    if (sw->spos < band_width) {
        for (i = 0; i < band_width - sw->spos; ++i) sw->F[i] = score_inf;
    }

    sw_align_sub(sw, 0, qpos);


    /* initialize the score matrix around the seed */
    int j_end;
    for (i = qpos; i < qpos + seedlen; ++i) {
        j_end = (i + 2) * colsize;
        for (j = (i + 1) * colsize; j < j_end; ++j) {
            sw->F[j] = score_inf;
        }

        assert(sw->query[i] == sw->subject[sw->spos + i + band_width]);

        j = (i + 1) * colsize + band_width;
        sw->F[j] = sw->F[j - colsize] + score_match;
    }


    /* align from the seed */
    sw_align_sub(sw, qpos + seedlen, qlen - 1);


    /* return the maximum score */
    sw->last_spos = 0;
    score_t s = score_inf;
    j_end = (qlen + 1) * colsize;

    if (spos + (qlen - qpos) + band_width >= sw->n) {
        j_end -= (spos + (qlen - qpos) + band_width) - sw->n;
    }

    for (j = qlen * colsize; j < j_end; ++j) {
        if (sw->F[j] < s) {
            /* TODO: check this shit */
            sw->last_spos = j - qlen * colsize;;
            s = sw->F[j];
        }
    }

    return s;
}


void sw_trace(sw_t* sw, sw_alignment_t* aln)
{
    edit_op_t op;
    int i, j, i_next, j_next;
    score_t s;

    aln->len = 0;

    for (i = sw->qlen, j = sw->last_spos; i > 0; i = i_next, j = j_next) {
        /* match / mismatch */
        i_next = i - 1;
        j_next = j;
        s  = sw->F[(i - 1) * colsize + j];
        if (sw->subject[sw->spos + (i - 1) + j] == sw->query[i - 1]) {
            op = EDIT_MATCH;
        }
        else {
            op = EDIT_MISMATCH;
        }

        /* s-gap */
        if (j < colsize - 1 &&
            sw->F[(i - 1) * colsize + j + 1] < s)
        {
            op = EDIT_S_GAP;
            s = sw->F[(i - 1) * colsize + j + 1];
            i_next = i - 1;
            j_next = j + 1;
        }

        /* q-gap */
        if (j > 0 &&
            sw->F[i * colsize + j - 1] < s)
        {
            op = EDIT_Q_GAP;
            s = sw->F[i * colsize + j - 1];
            i_next = i;
            j_next = j - 1;
        }


        /* make space for the alignment when needed */
        if (aln->size < aln->len + 1) {
            if (aln->size == 0) aln->size = 16;
            else aln->size *= 2;
            aln->ops = realloc_or_die(aln->ops, aln->size * sizeof(edit_op_t));
        }

        aln->ops[aln->len++] = op;
    }

    aln->spos = sw->spos + j - band_width;

    /* reverse the order of edit operations */
    i = 0;
    j = aln->len - 1;
    while (i < j) {
        op = aln->ops[i];
        aln->ops[i] = aln->ops[j];
        aln->ops[j] = op;
        ++i;
        --j;
    }
}





