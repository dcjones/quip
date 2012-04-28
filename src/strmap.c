
#include "strmap.h"
#include "misc.h"
#include <string.h>
#include <stdlib.h>

typedef struct strmap_value_t_
{
    str_t*   val;
    uint32_t idx;
} strmap_value_t;


struct strmap_t_
{
    strmap_value_t* xs;
    size_t m;           /* occupied positions */
    size_t n;           /* table size */
    size_t pn;          /* n = primes[pn] */
    size_t max_m;       /* max hashed items before size is doubled */
    size_t min_m;       /* min hashed items before size is halved */
};


/* Prime numbers that a near powers of two, suitable for hash table sizes, when
 * using quadratic probing. */
#define NUM_PRIMES 28
static const uint32_t primes[NUM_PRIMES] = {
           53U,         97U,        193U,        389U,    
          769U,       1543U,       3079U,       6151U,  
        12289U,      24593U,      49157U,      98317U,  
       196613U,     393241U,     786433U,    1572869U, 
      3145739U,    6291469U,   12582917U,   25165843U,
     50331653U,  100663319U,  201326611U,  402653189U,
    805306457U, 1610612741U, 3221225473U, 4294967291U };

static const double MAX_LOAD  = 0.7;
static const double MIN_LOAD  = 0.1;


/* simple quadratic probing */
static uint32_t probe(uint32_t h, uint32_t i)
{
    static const uint32_t c1 = 2;
    static const uint32_t c2 = 2;

    return h + i/c1 + (i*i)/c2;
}


strmap_t* strmap_alloc()
{
    strmap_t* M = malloc_or_die(sizeof(strmap_t));
    M->pn = 0;
    M->n  = primes[M->pn];
    M->xs = malloc_or_die(M->n * sizeof(strmap_value_t));
    memset(M->xs, 0, M->n * sizeof(strmap_value_t));
    M->m = 0;
    M->max_m = (size_t) ((double) M->n * MAX_LOAD);
    M->min_m = (size_t) ((double) M->n * MIN_LOAD);

    return M;
}


void strmap_free(strmap_t* M)
{
    if (M) {
        size_t i;
        for (i = 0; i < M->n; ++i) {
            str_free(M->xs[i].val);
            free(M->xs[i].val);
        }
        free(M->xs);
        free(M);
    }
}


size_t strmap_size(const strmap_t* M)
{
    return M->m;
}


static void ins_shallow(strmap_t* M, strmap_value_t* key)
{
    uint32_t h = murmurhash3(key->val->s, key->val->n);
    uint32_t probe_num = 1;
    uint32_t k = h % M->n;

    while (true) {
        if (M->xs[k].val == NULL) {
            break;
        }

        k = probe(h, ++probe_num) % M->n;
    }

    M->xs[k].val = key->val;
    M->xs[k].idx = key->idx;
}


static void strmap_expand(strmap_t* M)
{
    size_t new_n = primes[M->pn + 1];
    strmap_value_t* new_xs = malloc_or_die(new_n * sizeof(strmap_value_t));
    memset(new_xs, 0, new_n * sizeof(strmap_value_t));

    /* rehash */
    size_t i;
    for (i = 0; i < M->n; ++i) {
        if (M->xs[i].val) {
            ins_shallow(M, &M->xs[i]);
        }
    }

    free(M->xs);
    M->xs = new_xs;
    M->pn++;
    M->n = new_n;
    M->min_m = (size_t) ((double) M->n * MIN_LOAD);
    M->max_m = (size_t) ((double) M->n * MAX_LOAD);
}


uint32_t strmap_get(strmap_t* M, const str_t* key)
{
    if (M->m + 1 >= M->max_m) strmap_expand(M);

    uint32_t h = murmurhash3(key->s, key->n);
    uint32_t probe_num = 1;
    uint32_t k = h % M->n;

    while (true) {
        if (M->xs[k].val == NULL) {
            break;
        }
        else if (M->xs[k].val->n == key->n &&
                 strcmp((char*) M->xs[k].val->s, (char*) key->s) == 0) {
            return M->xs[k].idx;
        }

        k = probe(h, ++probe_num) % M->n;
    }

    M->xs[k].val = malloc_or_die(sizeof(str_t));
    str_copy(M->xs[k].val, key);
    M->xs[k].idx = M->m++;

    return M->xs[k].idx;
}

