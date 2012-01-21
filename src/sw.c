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


typedef uint16_t score_t;


struct sw_t_
{
    uint8_t* subject;
    uint8_t* query;

    int n; /* subject length */
    int m; /* query length */

    int qlen;
    int last_spos;

    /* minimum cost alignment matricies */

    /* alignments ending in some number of matches */
    score_t* M;

    /* alignments ending in some number of subject gaps */
    score_t* S;

    /* alignments ending in some number of query gaps */
    score_t* Q;

    /* overall minimum cost alignment */
    score_t* F;
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

static const score_t score_inf        = UINT16_MAX / 2;
static const score_t score_q_gap_open = 1;
static const score_t score_q_gap_ext  = 1;
static const score_t score_s_gap_open = 2;
static const score_t score_s_gap_ext  = 2;
static const score_t score_match_open = 1;
static const score_t score_match_ext  = 1;
static const score_t score_mismatch   = 2;




static inline score_t score_min(score_t a, score_t b)
{
    return a < b ? a : b;
}

static inline score_t score_min4(score_t a, score_t b, score_t c, score_t d)
{
    return score_min(score_min(a,b), score_min(c,d));
}


sw_t* sw_alloc(const twobit_t* subject)
{
    sw_t* sw = malloc_or_die(sizeof(sw_t));

    sw->n = (int) twobit_len(subject);

    sw->subject = malloc_or_die(sw->n * sizeof(uint8_t));

    int i;
    for (i = 0; i < sw->n; ++i) sw->subject[i] = twobit_get(subject, i);

    /* We don't know the query length in advance, so these are all allocated
     * and/or expanded when needed. */
    sw->m = 0;
    sw->query = NULL;
    sw->M = NULL;
    sw->S = NULL;
    sw->Q = NULL;
    sw->F = NULL;

    return sw;
}


void sw_free(sw_t* sw)
{
    free(sw->subject);
    free(sw->query);
    free(sw->M);
    free(sw->S);
    free(sw->Q);
    free(sw->F);
    free(sw);
}




static void sw_align_sub(sw_t* sw, int su, int sv, int qu, int qv)
{
    /* for convenience */
    const int k = sw->m + 1;

    int i, j, idx;
    for (i = su + 1; i <= sv + 1; ++i) {
        idx = i * k + qu + 1;

        for (j = qu + 1; j <= qv + 1; ++j, ++idx) {
            if (sw->subject[i - 1] == sw->query[j - 1]) {
                sw->M[idx] = score_min(sw->M[(i - 1) * k + (j - 1)] + score_match_ext,
                                       sw->F[(i - 1) * k + (j - 1)] + score_match_open);
            }
            else {
                sw->M[idx] = score_inf;
            }

            sw->Q[idx] = score_min(sw->Q[(i - 1) * k + j] + score_q_gap_ext,
                                   sw->F[(i - 1) * k + j] + score_q_gap_open);

            sw->S[idx] = score_min(sw->S[i * k + (j - 1)] + score_s_gap_ext,
                                   sw->F[i * k + (j - 1)] + score_s_gap_open);

            sw->F[idx] = score_min4(sw->F[(i - 1) * k + (j - 1)] + score_mismatch,
                                    sw->M[idx], sw->Q[idx], sw->S[idx]);
        }
    }
}


static void sw_ensure_query_len(sw_t* sw, int m)
{
    if (sw->m < m)  {
        sw->m = m;
        sw->query = realloc_or_die(sw->query, sw->m * sizeof(uint8_t));
        sw->M = realloc_or_die(sw->M, (sw->n + 1) * (sw->m + 1) * sizeof(score_t));
        sw->S = realloc_or_die(sw->S, (sw->n + 1) * (sw->m + 1) * sizeof(score_t));
        sw->Q = realloc_or_die(sw->Q, (sw->n + 1) * (sw->m + 1) * sizeof(score_t));
        sw->F = realloc_or_die(sw->F, (sw->n + 1) * (sw->m + 1) * sizeof(score_t));
    }
}




int sw_seeded_align(sw_t* sw, const twobit_t* query,
                    int spos, int qpos, int seedlen)
{
    int i, qlen = (int) twobit_len(query);
    sw->qlen = qlen;

    if (spos + seedlen + (qlen - qpos - seedlen + 1) >= sw->n) return -1;

    sw_ensure_query_len(sw, qlen);

    for (i = 0; i <= qlen; ++i) {
        sw->query[i] = twobit_get(query, i);
    }


    /* Align up to the seed.
     * =====================
     * */

    /* We limit of the amount of subject sequence consumed in the alignment,
     * to make the run time O(m^2), rather than O(mn) */
    int s0 = spos - 2 * qpos;
    if (s0 < 0) s0 = 0;

    int s1 = spos + seedlen + 2 * (qlen - qpos - seedlen + 1);
    if (s1 >= sw->n) s1 = sw->n - 1;


    /* initialize column 0 */
    int idx;
    for (i = s0; i <= spos + 1; ++i) {
        idx = i * (sw->m + 1);
        sw->M[idx] = score_inf;
        sw->S[idx] = score_inf;
        sw->Q[idx] = score_inf;
        sw->F[idx] = 0;
    }

    /* initialize row 0 */
    for (i = 1; i < qpos + 1; ++i) {
        idx = s0 * (sw->m + 1) + i;
        sw->M[idx] = score_inf;
        sw->S[idx] = score_inf;
        sw->Q[idx] = score_q_gap_open + (i - 1) * score_q_gap_ext;
        sw->F[idx] = sw->Q[idx];
    }

    sw_align_sub(sw, s0, spos, 0, qpos);


    /* Set the score matrix around the seed */
    int idx2;
    for (i = 0; i < seedlen; ++i) {
        assert(sw->query[qpos + i] == sw->subject[spos + i]);

        idx = (spos + i) * (sw->m + 1) + (qpos + i);
        idx2 = (spos + 1 + i) * (sw->m + 1) + (qpos + 1 + i);
        sw->M[idx2] = sw->F[idx2] = sw->F[idx] + score_match_ext;

        idx2 = (spos + 1 + i + 1) * (sw->m + 1) + (qpos + 1 + i);
        sw->F[idx2] = score_inf;

        idx2 = (spos + 1 + i) * (sw->m + 1) + (qpos + 1 + i + 1);
        sw->F[idx2] = score_inf;
    }

    idx = (spos + i) * (sw->m + 1) + (qpos + i);
    idx2 = (spos + 1 + i) * (sw->m + 1) + (qpos + 1 + i);
    sw->M[idx2] = sw->F[idx2] = sw->F[idx] + score_match_ext;
    assert(sw->F[idx2] < score_inf);


    /* Align from the seed.
     * ====================
     */

    /* initialize column 0 */
    for (i = spos + seedlen + 1; i <= s1 + 1; ++i) {
        idx = i * (sw->m + 1) + qpos + seedlen;
        sw->M[idx] = score_inf;
        sw->S[idx] = score_inf;
        sw->Q[idx] = sw->F[idx2] + score_q_gap_open + (i - (spos + seedlen)) * score_q_gap_ext;
        sw->F[idx] = sw->Q[idx];
    }

    /* initialize row 0 */
    for (i = qpos + seedlen + 1; i < qlen + 1; ++i) {
        idx = (spos + seedlen) * (sw->m + 1) + i;
        sw->M[idx] = score_inf;
        sw->S[idx] = sw->F[idx2] + score_s_gap_open + (i - (qpos + seedlen + 1)) * score_s_gap_ext;
        sw->Q[idx] = score_inf;
        sw->F[idx] = sw->S[idx];
    }

    sw_align_sub(sw, spos + seedlen, s1, qpos + seedlen, qlen - 1);

    score_t s = score_inf;
    for (i = spos + seedlen; i <= s1; ++i) {
        if (sw->F[i * (sw->m + 1) + qlen] < s) {
            s = sw->F[i * (sw->m + 1) + qlen];
            sw->last_spos = i;
        }
        s = score_min(s, sw->F[i * (sw->m + 1) + qlen]);
    }

    return s;
}


void sw_trace(sw_t* sw, sw_alignment_t* aln)
{
    aln->len = 0;
    int i, j;
    int i_next, j_next;
    edit_op_t op;

    i = sw->last_spos;
    j = sw->qlen;

    score_t s = sw->F[i * (sw->m + 1) + j];

    while (j > 0) {
        if (aln->size < aln->len + 1) {
            if (aln->size == 0) aln->size = 16;
            else aln->size *= 2;
            aln->ops = realloc_or_die(aln->ops, aln->size * sizeof(edit_op_t));
        }

        s = sw->F[i * (sw->m + 1) + j - 1];
        i_next = i;
        j_next = j - 1;
        op = EDIT_S_GAP;

        if (i > 0) {
            if (sw->F[(i - 1) * (sw->m + 1) + j] <= s) {
                s = sw->F[(i - 1) * (sw->m + 1) + j];
                i_next = i - 1;
                j_next = j;
                op = EDIT_Q_GAP;
            }

            if (sw->F[(i - 1) * (sw->m + 1) + (j - 1)] <= s) {
                s = sw->F[(i - 1) * (sw->m + 1) + (j - 1)];
                i_next = i - 1;
                j_next = j - 1;

                /* match or mismatch */
                op = sw->subject[i - 1] == sw->query[j - 1] ? EDIT_MATCH : EDIT_MISMATCH;
            }
        }

        aln->ops[aln->len++] = op;

        i = i_next;
        j = j_next;
    }

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

    aln->spos = i;
}





