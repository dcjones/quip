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


typedef enum edit_op_t_
{
    EDIT_MATCH,
    EDIT_MISMATCH,
    EDIT_Q_GAP,
    EDIT_S_GAP
} edit_op_t;


struct sw_t_
{
    uint8_t* subject;
    uint8_t* query;
    int size;

    /* matrix columns */
    int* ws;
    int* ys;
    int* wxyzs;

    /* path */
    edit_op_t* ps;
    int path_len;
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


    sw->ps = malloc_or_die(sw->size * sizeof(edit_op_t));
    sw->path_len = 0;

    return sw;
}


void sw_free(sw_t* sw)
{
    free(sw->subject);
    free(sw->query);
    free(sw->ws);
    free(sw->ys);
    free(sw->wxyzs);
    free(sw->ps);
    free(sw);
}



/* Forward alignment restricted to [su, sv] in the subject sequence and [qu, qv]
 * in the query sequence. If 'anchor' is true, the alignment must begin with a
 * match or mismatch. */
static int sw_align_foreward(sw_t* sw, int su, int sv, int qu, int qv, bool anchor)
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
        ws[su] = imin(wxyzs[su] + score_s_gap_open, ws[su] + score_s_gap_ext);

        /* query gap */
        x = score_pos_inf;

        /* mismatch */
        z = score_mismatch;

        /* match */
        if (j == qu) {
            if (sw->subject[su] == sw->query[qu]) {
                ys[su] = score_match_open;
            }
            else {
                ys[su] = score_pos_inf;
            }
        }
        else {
            ys[su] = score_pos_inf;
        }

        wxyz0 = wxyzs[su];

        if (j == qu && anchor) {
            wxyzs[su] = imin(ys[su], z);
        }
        else {
            wxyzs[su] = imin4(ws[su], x, ys[su], z);
        }

        y0 = score_pos_inf;
        wxyz_min = imin(wxyz_min, wxyzs[su]);


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



static int sw_align_backward(sw_t* sw, int su, int sv, int qu, int qv, bool anchor)
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


    for (j = qv; j >= qu; --j) {

        wxyz_min = score_pos_inf;

        /* last row special case ? */

        /* subject gap */
        ws[sv] = imin(wxyzs[sv] + score_s_gap_open, ws[sv] + score_s_gap_ext);

        /* query gap */
        x = score_pos_inf;

        /* mismatch */
        z = score_mismatch;

        /* match */
        if (j == qv) {
            if (sw->subject[sv] == sw->query[qv]) {
                ys[sv] = score_match_open;
            }
            else {
                ys[sv] = score_pos_inf;
            }
        }
        else {
            ys[sv] = score_pos_inf;
        }

        wxyz0 = wxyzs[sv];

        if (j == qv && anchor) {
            wxyzs[sv] = imin(ys[sv], z);
        }
        else {
            wxyzs[sv] = imin4(ws[sv], x, ys[sv], z);
        }

        y0 = score_pos_inf;
        wxyz_min = imin(wxyz_min, wxyzs[sv]);


        for (i = sv - 1; i >= su; --i) {
            /* subject gap */
            ws[i] = imin(wxyzs[i] + score_s_gap_open, ws[i] + score_s_gap_ext);

            /* query gap */
            x = imin(wxyzs[i + 1] + score_q_gap_open, x + score_q_gap_ext);

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


static void sw_align_path(sw_t* sw,
                          int su, int sv, int qu, int qv,
                          bool anchor00, bool anchor11)
{
    int q;

    /* TODO:
     * Can we actually divide and conqueror with sv - su == 1,
     * or with qv - qu == 1?
     */

    if (sv - su == 0) {
        // TODO
        return;
    }

    if (qv - qu == 0) { 
        // TODO
        return;
    }

    sw_align_foreward(sw, su, sv, qu, qu + (qv - qu) / 2, anchor00);

    // TODO: do something with wxyzs

    sw_align_backward(sw, su, sv, qu + (qv - qu) / 2 + 1, qv, anchor11);

    // TODO: find 'q'

    // TODO: push q to sw->ps

    sw_align_path(sw, su, q, qu, qu + (qv - qu) / 2, anchor00, true);
    sw_align_path(sw, q, sv, qu + (qv - qu) / 2  + 1, qv, true, anchor11);

    // TODO: rearrange the path, perhaps ?
}



int sw_align(sw_t* sw, const twobit_t* query,
             int spos, int qpos, int seedlen)
{
    size_t i, n = twobit_len(query);
    for (i = 0; i < n; ++i) {
        sw->query[i] = twobit_get(query, i);
    }

    /* we bound the extent of the alignment within the subject sequence, making
     * the algorithm O(m^2) rather than O(mn)
     */

    /* least allowable subject position */
    int s0 = spos - 2 * qpos;
    if (s0 < 0) s0 = 0;

    /* greatest allowable subject position */
    int s1 = spos + seedlen + 2 * ((int) n - qpos + seedlen);
    if (s1 >= sw->size) s1 = sw->size - 1;

    /* align sequence left of the seed */
    sw_align_path(sw, s0, spos, 0, qpos, false, true);

    /* align to the right of the seed */
    sw_align_path(sw, spos + seedlen - 1, s1, qpos + seedlen - 1, (int) n - 1, true, false);

    /* TODO score */

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

    sw_align_foreward(sw, 0, sw->size - 1, 0, n - 1, false);

    return 0;
}



