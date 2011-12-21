
#include "qual.h"
#include "misc.h"
#include <stdint.h>
#include <stdlib.h>
#include <assert.h>


/* Every quality encoding scheme uses ASCII charocters in [33, 126] */
const char   qual_first = 33;
const char   qual_last  = 126;
const size_t qual_size  = 94;


/* We expand arrays as we see longer reads. This is the expansion regime, which
 * is based on common read lengths*/
const size_t ks[] = {36, 50, 55, 75, 100, 150};
const size_t ks_size = sizeof(ks);


struct qualmodel_t_
{
    /* Score given position and nucleotide.
     * Indexed as: position -> (score x nucleotide)
     */
    uint32_t** prS;

    /* Run-length , given position and score. 
     * Indexed as: position -> (score x run-length)
     */
    uint32_t** prK;

    /* Maximum read length .*/
    size_t k;
};

static void setallones(uint32_t* xs, size_t n)
{
    size_t i;
    for (i = 0; i < n; ++i) xs[i] = 1;
}


qualmodel_t* qualmodel_alloc()
{
    qualmodel_t* M = malloc_or_die(sizeof(qualmodel_t));

    M->k = ks[0];

    size_t i;

    M->prS = malloc_or_die(M->k * sizeof(uint32_t*));
    M->prK = malloc_or_die(M->k * sizeof(uint32_t*));
    for (i = 0; i < M->k; ++i) {
        M->prS[i] = malloc_or_die(qual_size * 5 * sizeof(uint32_t));
        setallones(M->prS[i], qual_size * 5);

        M->prK[i] = malloc_or_die(qual_size * (M->k - i) * sizeof(uint32_t));
        setallones(M->prK[i], qual_size * (M->k - i));
    }

    return M;
}


void qualmodel_free(qualmodel_t* M)
{
    size_t i;
    for (i = 0; i < M->k; ++i) {
        free(M->prS[i]);
        free(M->prK[i]);
    }
    free(M->prS);
    free(M->prK);
}



/* Expand the model to handle read length of at least min_k. */
static void quadmodel_expand(qualmodel_t* M, size_t min_k)
{
    size_t new_k;
    
    /* use a predetermined size if there is one */
    size_t i = 0;
    do {
        new_k = ks[i];
        ++i;
    } while(i < ks_size && new_k < min_k);

    while (new_k < min_k) new_k = min_k;

    assert(new_k > M->k);


    M->prS = realloc_or_die(M->prS, new_k * sizeof(uint32_t*));
    for (i = M->k; i < new_k; ++i) {
        M->prS[i] = malloc_or_die(qual_size * 5 * sizeof(uint32_t));
        setallones(M->prS[i], qual_size * 5);
    }

    for (i = 0; i < M->k; ++i) {
        M->prK[i] = realloc_or_die(M->prK[i], qual_size * (new_k - i) * sizeof(uint32_t));
        setallones(M->prK[i] + qual_size * (M->k - i), qual_size * (new_k - M->k));
    }

    M->prK = realloc_or_die(M->prK, new_k * sizeof(uint32_t*));
    for (i = M->k; i < new_k; ++i) {
        M->prK[i] = malloc_or_die(qual_size * (new_k - i) * sizeof(uint32_t));
        setallones(M->prK[i], qual_size * (new_k - i));
    }

    M->k = new_k;
}


unsigned char ascii_to_nucnum(char x)
{
    switch (x) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default:  return 4;
    }
}


/* TODO: this function should just work on one sequence */
void qualmodel_update(qualmodel_t* M, const seq_t* x)
{
    unsigned char q; // quality score
    unsigned char u; // nucleotide
    size_t ki, kj;   // position
    size_t l;        // run-length

    if (x->qual.n > M->k) quadmodel_expand(M, x->qual.n);

    u = ascii_to_nucnum(x->seq.s[0]);
    q = (unsigned char) (x->qual.s[0] - qual_first);
    M->prS[0][u * qual_size + q] += 1;

    ki = 0;
    kj = 1;
    while (kj < x->qual.n) {
        u = ascii_to_nucnum(x->seq.s[kj]);
        q = (unsigned char) (x->qual.s[kj] - qual_first);
        M->prS[kj][u * qual_size + q] += 1;

        if (x->qual.s[kj] != x->qual.s[ki]) {
            l = kj - ki;
            q = (unsigned char) (x->qual.s[ki] - qual_first);
            M->prK[ki][(l - 1) * qual_size + q] += 1;

            ki = kj;
        }

        ++kj;
    }

    l = kj - ki;
    q = (unsigned char) (x->qual.s[ki] - qual_first);
    M->prK[ki][(l - 1) * qual_size + q] += 1;
}

