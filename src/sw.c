/*
 * This file is part of quip.
 *
 * Copyright (c) 2011 by Daniel C. Jones <dcjones@cs.washington.edu>
 *
 * This is a very simple, 'fast enough' implementation of the Smith-Waterman
 * algorithm specifically for short nucleotide sequences, working in O(mn) time
 * and O(m) space, implemented according to the original Gotoh paper and
 * Phil Green's implementation in cross_match.
 *
 */


#include "sw.h"
#include "misc.h"
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <limits.h>


struct sw_t_
{
    uint8_t* subject;
    uint8_t* query;
    int size;

    /* matrix columns, used internally */
    int* ws;
    int* ys;
    int* wxyzs;
};


/* The alignment algorithm we use has a very add scoring system.  We assign
 * scores that correspond to the number of bits needed to represent the query
 * sequence, assuming the subject sequence is known.
 *
 * This way, the alignment reported is the alignment that minimizes
 * representation space, rather than simularity, ro edit-distance.
 *
 * Below are estimates of the number of bits needed, assuming a simple run
 * length encoding scheme is used.
 */

/* TODO: a more realisitic scoring system */
static const int score_neg_inf    = INT_MIN / 2;
static const int score_pos_inf    = INT_MAX / 2;
static const int score_q_gap_open = 16;
static const int score_q_gap_ext  = 0;
static const int score_s_gap_open = 8;
static const int score_s_gap_ext  = 8;
static const int score_match_open = 16;
static const int score_match_ext  = 0;
static const int score_mismatch   = 8;




static inline int imin(int a, int b)
{
    return a < b ? a : b;
}

static inline int imin4(int a, int b, int c, int d)
{
    return imin(imin(a,b), imin(c,d));
}


sw_t* sw_alloc(const twobit_t* subject)
{
    sw_t* sw = malloc_or_die(sizeof(sw_t));

    sw->size = (int) twobit_len(subject);

    sw->subject = malloc_or_die(sw->size * sizeof(uint8_t));

    /* In this particular application, query sequences can never be longer that
     * the subject sequence, so we reserve */
    sw->query = malloc_or_die(sw->size * sizeof(uint8_t));

    int i;
    for (i = 0; i < sw->size; ++i) sw->subject[i] = twobit_get(subject, i);

    sw->ws    = malloc_or_die(sw->size * sizeof(int));
    sw->ys    = malloc_or_die(sw->size * sizeof(int));
    sw->wxyzs = malloc_or_die(sw->size * sizeof(int));

    return sw;
}


void sw_free(sw_t* sw)
{
    free(sw->subject);
    free(sw->query);
    free(sw->ws);
    free(sw->ys);
    free(sw->wxyzs);
    free(sw);
}



/* Forward alignment restricted to [su, sv] in the subject sequence and [qu, qv]
 * in the query sequence.  */
static int sw_align_foreward(sw_t* sw, int su, int sv, int qu, int qv)
{
    int i, j;

    /* for convenience */
    int* ws = sw->ws;
    int x;
    int* ys = sw->ys;
    int z;
    int* wxyzs = sw->wxyzs;

    int y0, y1, wxyz0, wxyz_min;


    /* initialize ws to prohibit openining with a subject gap extension */
    for (i = su; i <= sv; ++i) ws[i] = score_pos_inf;

    /* initialize ys to prohibit opening with a match extension */
    for (i = su; i < sv; ++i) ys[i] = score_pos_inf;

    /* initialize wxyzs to allow starting from any position in the subject
     * sequence */
    memset(wxyzs + su, 0, (sv - su + 1) * sizeof(int));


    for (j = qu; j <= qv; ++j) {

        wxyz_min = score_pos_inf;

        /* first row special case */

        /* subject gap */
        wxyzs[su] = ws[su] = imin(wxyzs[su] + score_s_gap_open, ws[su] + score_s_gap_ext);

        /* query gap */
        x = score_pos_inf;

        /* match */
        if (j == qu) {
            if (sw->subject[su] == sw->query[qu]) {
                ys[su] = score_match_open;
            }
            else {
                ys[su] = score_pos_inf;
                z = score_mismatch;
            }
        }
        else {
            ys[su] = score_pos_inf;
            z      = score_pos_inf;
        }

        wxyz0 = wxyzs[su];
        wxyzs[su] = imin4(ws[su], x, ys[su], z);

        y0    = score_pos_inf;

        wxyz_min = imin(wxyz_min, wxyzs[su]);


        /* prohibit leading with a query gap extension */
        x = score_pos_inf;

        for (i = su + 1; i <= sv; ++i) {
            /* subject gap */
            ws[i] = imin(wxyzs[i] + score_s_gap_open, ws[i] + score_s_gap_ext);

            /* query gap */
            x = imin(wxyzs[i - 1] + score_q_gap_open, x + score_q_gap_ext);

            /* match */
            y1 = ys[i];
            if (sw->subject[i] == sw->query[j]) {
                ys[i] = imin(wxyz0 + score_match_open, y0 + score_match_ext);
            }
            else ys[i] = score_pos_inf;
            y0 = y1;

            /* mismatch */
            z = wxyz0 + score_mismatch;

            wxyz0 = wxyzs[i];
            wxyzs[i] = imin4(ws[i], x, ys[i], z);
            wxyz_min = imin(wxyz_min, wxyzs[i]);
        }
    }

    return wxyz_min;
}



static int sw_align_backward(sw_t* sw, int su, int sv, int qu, int qv)
{
    /* TODO */

    return 0;
}



int sw_align(sw_t* sw, const twobit_t* query,
             int spos, int qpos, int seedlen)
{
    size_t i, n = twobit_len(query);
    for (i = 0; i < n; ++i) {
        sw->query[i] = twobit_get(query, i);
    }

    // XXX: just for testing
    sw_align_foreward(sw, 0, sw->size - 1, 0, n - 1);


    /* TODO: align up to (spos, qpos) */

    /* TODO: align from (spos + seedlen, qpos + seedlen)
     * to (x, n), where x is anything.
     */

    return 0;
}




sw_t* sw_alloc_char(const char* subject)
{
    sw_t* sw = malloc_or_die(sizeof(sw_t));

    sw->size = (int) strlen(subject);

    sw->subject = malloc_or_die(sw->size * sizeof(uint8_t));

    /* In this particular application, query sequences can never be longer that
     * the subject sequence, so we reserve */
    sw->query = malloc_or_die(sw->size * sizeof(uint8_t));

    memcpy(sw->subject, subject, sw->size);

    sw->ws    = malloc_or_die(sw->size * sizeof(int));
    sw->ys    = malloc_or_die(sw->size * sizeof(int));
    sw->wxyzs = malloc_or_die(sw->size * sizeof(int));

    return sw;

}


int sw_align_char(sw_t* sw, const char* query)
{
    size_t n = strlen(query);
    memcpy(sw->query, query, n);

    sw_align_foreward(sw, 0, sw->size - 1, 0, n - 1);
}



