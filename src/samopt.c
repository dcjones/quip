

#include "quip.h"
#include "kmer.h"
#include "misc.h"
#include <string.h>


struct samopt_table_t_
{
    samopt_t* xs;
    size_t m;           /* occupied positions */
    size_t n;           /* table size */
    size_t pn;          /* n = primes[pn] */
    size_t max_m;       /* max hashed items before size is doubled */

};


static uint32_t samopt_key_hash(const unsigned char key[2])
{
    /* TODO: a more appropriate hash function */
    return kmer_hash(*(kmer_t*) key);
}


static bool samopt_table_empty(const samopt_t* opt)
{
    return opt->key[0] == '\0' && opt->key[1] == '\0';
}


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

/* simple quadratic probing */
static uint32_t probe(uint32_t h, uint32_t i)
{
    static const uint32_t c1 = 2;
    static const uint32_t c2 = 2;

    return h + i/c1 + (i*i)/c2;
}


samopt_table_t* samopt_table_alloc()
{
    samopt_table_t* M = malloc_or_die(sizeof(samopt_table_t));
    M->pn = 0;
    M->n  = primes[M->pn];
    M->xs = malloc_or_die(M->n * sizeof(samopt_t));
    memset(M->xs, 0, M->n * sizeof(samopt_t));
    M->m = 0;
    M->max_m = (size_t) ((double) M->n * MAX_LOAD);

    return M;
}


void samopt_table_free(samopt_table_t* M)
{
    if (M) {
        size_t i;
        for (i = 0; i < M->n; ++i) {
            str_free(M->xs[i].data);
        }
        free(M->xs);
        free(M);
    }
}


size_t samopt_table_size(const samopt_table_t* M)
{
    return M->m;
}


static void samopt_table_expand(samopt_table_t* M)
{
    size_t old_n = M->n;
    samopt_t* old_xs = M->xs;

    M->n = primes[++M->pn];
    M->xs = malloc_or_die(M->n * sizeof(samopt_t));
    memset(M->xs, 0, M->n * sizeof(samopt_t));

    /* rehash */
    size_t i;
    uint32_t h, probe_num, k;

    for (i = 0; i < old_n; ++i) {
        /* shallow insert */
        if (samopt_table_empty(&old_xs[i])) continue;

        h = samopt_key_hash(old_xs[i].key);
        probe_num = 1;
        k = h % M->n;

        while (true) {
            if (samopt_table_empty(&M->xs[k])) break;
            k = probe(h, ++probe_num) % M->n;
        }

        M->xs[k].key[0] = old_xs[i].key[0];
        M->xs[k].key[1] = old_xs[i].key[1];
        M->xs[k].type = old_xs[i].type;
        M->xs[k].data = old_xs[i].data;
    }

    free(old_xs);
    M->max_m = (size_t) ((double) M->n * MAX_LOAD);
}


static samopt_t* samopt_table_get_priv(const samopt_table_t* M, const unsigned char key[2])
{
    uint32_t h = samopt_key_hash(key);
    uint32_t probe_num = 1;
    uint32_t k = h & M->n;

    while (true) {
        if (samopt_table_empty(&M->xs[k]) ||
            (M->xs[k].key[0] == key[0] && M->xs[k].key[1] == key[1]))
        {
            return &M->xs[k];
        }

        k = probe(h, ++probe_num) & M->n;
    }

    return NULL;
}


samopt_t* samopt_table_get(samopt_table_t* M, const unsigned char key[2])
{
    if (M->m + 1 >= M->max_m) samopt_table_expand(M);

    samopt_t* opt = samopt_table_get_priv(M, key);

    if (samopt_table_empty(opt)) {
        opt->key[0] = key[0];
        opt->key[1] = key[1];
    }

    if (opt->data == NULL) {
        opt->data = malloc_or_die(sizeof(str_t));
        str_init(opt->data);
    }

    return opt;
}

#if 0
void samopt_table_put(samopt_table_t* M, unsigned char key[2], unsigned char type, const str_t* data)
{
    if (M->m + 1 >= M->max_m) samopt_table_expand(M);

    samopt_t* opt = samopt_table_get_priv(M, key);
    opt->key[0] = key[0]; opt->key[1] = key[1];
    opt->type = type;
    if (opt->data == NULL) {
        opt->data = malloc_or_die(sizeof(str_t));
        str_init(opt->data);
    }
    str_copy(opt->data, data);
}
#endif


void samopt_table_clear(samopt_table_t* M)
{
    size_t i;
    for (i = 0; i < M->n; ++i) {
        if (M->xs[i].data) M->xs[i].data->n = 0;
    }

    M->m = 0;
}


