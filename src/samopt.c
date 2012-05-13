

#include "samoptenc.h"
#include "kmer.h"
#include "misc.h"
#include "ac.h"
#include "dist.h"
#include <string.h>
#include <assert.h>


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
    kmer_t k = ((kmer_t) key[1] << 8) | (kmer_t) key[0];
    return kmer_hash(k);
}


static bool samopt_table_empty(const samopt_t* opt)
{
    return opt->key[0] == '\0' && opt->key[1] == '\0';
}


/* Prime numbers that a near powers of two, suitable for hash table sizes, when
 * using quadratic probing. */
#define NUM_PRIMES 30
static const uint32_t primes[NUM_PRIMES] = {
                                     11U,         23U,
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


size_t samopt_table_bytes(const samopt_table_t* M)
{
    size_t i;
    size_t bytes = 0;
    for (i = 0; i < M->n; ++i) {
        if (samopt_table_empty(&M->xs[i])) continue;
        else if (M->xs[i].data && M->xs[i].data->n > 0) {
            bytes += 5 + M->xs[i].data->n;
        }
    }

    /* don't count a trailing tab */
    bytes -= 1;

    return bytes;
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
    uint32_t k = h % M->n;

    while (true) {
        if (samopt_table_empty(&M->xs[k]) ||
            (M->xs[k].key[0] == key[0] && M->xs[k].key[1] == key[1]))
        {
            return &M->xs[k];
        }

        k = probe(h, ++probe_num) % M->n;
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
        M->m++;
    }

    if (opt->data == NULL) {
        opt->data = malloc_or_die(sizeof(str_t));
        str_init(opt->data);
    }

    return opt;
}


void samopt_table_clear(samopt_table_t* M)
{
    size_t i;
    for (i = 0; i < M->n; ++i) {
        if (M->xs[i].data) M->xs[i].data->n = 0;
    }
}

void samopt_table_copy(samopt_table_t* dest, const samopt_table_t* src)
{
    samopt_table_clear(dest);
    samopt_t* opt;
    size_t i;
    for (i = 0; i < src->n; ++i) {
        if (samopt_table_empty(&src->xs[i])) continue;

        opt = samopt_table_get(dest, src->xs[i].key);
        opt->type = src->xs[i].type;
        str_copy(opt->data, src->xs[i].data);
    }
}


void samopt_table_append_str(const samopt_table_t* M, str_t* dest)
{
    size_t i;
    for (i = 0; i < dest->n; ++i) {
        if (samopt_table_empty(&M->xs[i]) ||
            M->xs[i].data == NULL ||
            M->xs[i].data->n == 0) continue;


        str_reserve_extra(dest, 7 + M->xs[i].data->n);
        dest->n += snprintf((char*) dest->s + dest->n, 6, "\t%c%c:%c:",
                            M->xs[i].key[0], M->xs[i].key[1], M->xs[i].type);
        str_append(dest, M->xs[i].data);
    }
}


/* the size of keyspace: [A-Za-z][A-Za-z0-9] */
#define keys_N 3224


typedef enum samoptenc_flag_t_ {
    SAMOPTENC_FLAG_INSERT,
    SAMOPTENC_FLAG_START
} samoptenc_flag_t;

struct samoptenc_t_
{
    ac_t* ac;

    /* a table containing all tags seen so far */
    samopt_table_t* last;

    /* distribution leading insertion flags */
    dist2_t d_encflag;

    /* distribution over tag characters */
    dist128_t d_tag_char;

    /* types, conditioned on tag */
    dist16_t* d_type[keys_N];

    /* data, conditioned on tag */
    cond_dist128_t* d_data[keys_N];
};


static void samoptenc_init(samoptenc_t* E)
{
    E->last = samopt_table_alloc();

    dist2_init(&E->d_encflag);
    dist128_init(&E->d_tag_char);
    memset(E->d_type, 0, keys_N * sizeof(dist16_t*));
    memset(E->d_data, 0, keys_N * sizeof(dist128_t*));
}


samoptenc_t* samoptenc_alloc_encoder(quip_writer_t writer, void* writer_data)
{
    samoptenc_t* E = malloc_or_die(sizeof(samoptenc_t));
    E->ac = ac_alloc_encoder(writer, writer_data);
    samoptenc_init(E);

    return E;
}


samoptenc_t* samoptenc_alloc_decoder(quip_reader_t reader, void* reader_data)
{
    samoptenc_t* E = malloc_or_die(sizeof(samoptenc_t));
    E->ac = ac_alloc_decoder(reader, reader_data);
    samoptenc_init(E);

    return E;
}


void samoptenc_free(samoptenc_t* E)
{
    if (E) {
        ac_free(E->ac);
        samopt_table_free(E->last);

        size_t i;
        for (i = 0; i < keys_N; ++i) {
            free(E->d_type[i]);
            if (E->d_data[i]) {
                cond_dist128_free(E->d_data[i]);
                free(E->d_data[i]);
            }
        }

        free(E);
    }
}


/* Map a key in [A-Za-z][A-Za-z0-9] to a number in [0,key_N] */
static uint32_t samopt_key_num(const unsigned char key[2])
{
    uint32_t x = 0;

    if      ('A' <= key[0] && key[0] <= 'Z') x = (uint32_t) (key[0] - 'A');
    else if ('a' <= key[0] && key[0] <= 'z') x = 26 + (uint32_t) (key[0] - 'a');
    x *= 52;

    if      ('A' <= key[1] && key[1] <= 'Z') x += (uint32_t) (key[1] - 'A');
    else if ('a' <= key[1] && key[1] <= 'z') x += 26 + (uint32_t) (key[1] - 'a');
    else if ('0' <= key[1] && key[1] <= '9') x += 52 + (uint32_t) (key[1] - '0');

    return x;
}


/* Map a field type into [0, 16] */
static uint32_t samopt_type_num(unsigned char t)
{
    switch (t) {
        case 'A': return  0;
        case 'C': return  1;
        case 'c': return  2;
        case 'S': return  3;
        case 's': return  4;
        case 'I': return  5;
        case 'i': return  6;
        case 'f': return  7;
        case 'd': return  8;
        case 'Z': return  9;
        case 'H': return 10;
        case 'B': return 11;
        default:
            quip_error("Unknown SAM field type: %c", t);
            return 16;
    }
}

static unsigned char samopt_type_num_rev(uint32_t t)
{
    switch (t) {
        case  0: return 'A';
        case  1: return 'C';
        case  2: return 'c';
        case  3: return 'S';
        case  4: return 's';
        case  5: return 'I';
        case  6: return 'i';
        case  7: return 'f';
        case  8: return 'd';
        case  9: return 'Z';
        case 10: return 'H';
        case 11: return 'B';
        default:
            quip_error("Unknown SAM field type: %c", t);
            return '\0';
    }
}


void samoptenc_encode(samoptenc_t* E, const samopt_table_t* T)
{
    samopt_table_clear(E->last);

    uint32_t k;
    samopt_t* opt;
    size_t i, j;
    for (i = 0; i < T->n; ++i) {
        if (samopt_table_empty(&T->xs[i])) continue;
        opt = samopt_table_get_priv(E->last, T->xs[i].key);

        if (samopt_table_empty(opt)) {
            opt->key[0] = T->xs[i].key[0];
            opt->key[1] = T->xs[i].key[1];
            E->last->m++;

            if (opt->data == NULL) {
                opt->data = malloc_or_die(sizeof(str_t));
                str_init(opt->data);
            }

            k = samopt_key_num(opt->key);

            E->d_type[k] = malloc_or_die(sizeof(dist16_t));
            dist16_init(E->d_type[k]);

            E->d_data[k] = malloc_or_die(sizeof(cond_dist128_t));
            cond_dist128_init(E->d_data[k], 128);

            dist2_encode(E->ac, &E->d_encflag, SAMOPTENC_FLAG_INSERT);
            dist128_encode(E->ac, &E->d_tag_char, T->xs[i].key[0]);
            dist128_encode(E->ac, &E->d_tag_char, T->xs[i].key[1]);
        }

        opt->type= T->xs[i].type;
        str_copy(opt->data, T->xs[i].data);
    }

    dist2_encode(E->ac, &E->d_encflag, SAMOPTENC_FLAG_START);

    for (i = 0; i < E->last->n; ++i) {
        if (samopt_table_empty(&E->last->xs[i])) continue;

        opt = &E->last->xs[i];

        k = samopt_key_num(opt->key);
        dist16_encode(E->ac, E->d_type[k], samopt_type_num(opt->type));

        unsigned char c = 0;
        for (j = 0; j < opt->data->n; ++j) {
            cond_dist128_encode(E->ac, E->d_data[k], c, opt->data->s[j]);
            c = opt->data->s[j];
        }

        cond_dist128_encode(E->ac, E->d_data[k], c, '\0');
    }
}


void samoptenc_decode(samoptenc_t* E, samopt_table_t* T)
{
    samopt_table_clear(E->last);

    samopt_t* opt;
    uint32_t k;
    while (dist2_decode(E->ac, &E->d_encflag) == SAMOPTENC_FLAG_INSERT) {
        unsigned char key[2];
        key[0] = dist128_decode(E->ac, &E->d_tag_char);
        key[1] = dist128_decode(E->ac, &E->d_tag_char);

        opt = samopt_table_get_priv(E->last, key);
        assert(samopt_table_empty(opt));

        opt->key[0] = key[0];
        opt->key[1] = key[1];
        E->last->m++;

        if (opt->data == NULL) {
            opt->data = malloc_or_die(sizeof(str_t));
            str_init(opt->data);
        }

        k = samopt_key_num(opt->key);

        E->d_type[k] = malloc_or_die(sizeof(dist16_t));
        dist16_init(E->d_type[k]);

        E->d_data[k] = malloc_or_die(sizeof(cond_dist128_t));
        cond_dist128_init(E->d_data[k], 128);
    }

    size_t i;
    for (i = 0; i < E->last->n; ++i) {
        if (samopt_table_empty(&E->last->xs[i])) continue;

        opt = &E->last->xs[i];

        k = samopt_key_num(opt->key);
        opt->type = samopt_type_num_rev(dist16_decode(E->ac, E->d_type[k]));

        unsigned char c = 0;
        while (true) {
            if (opt->data->n + 2 >= opt->data->size) str_reserve_extra(opt->data, 16);
            c = opt->data->s[opt->data->n] = cond_dist128_decode(E->ac, E->d_data[k], c);

            if (c == '\0') break;
            else opt->data->n++;
        }

        opt->data->s[opt->data->n] = '\0';
    }

    samopt_table_copy(T, E->last);
}


size_t samoptenc_finish(samoptenc_t* E)
{
    return ac_finish_encoder(E->ac);
}


void samoptenc_flush(samoptenc_t* E)
{
    ac_flush_encoder(E->ac);
}


void samoptenc_start_decoder(samoptenc_t* E)
{
    ac_start_decoder(E->ac);
}


void samoptenc_reset_decoder(samoptenc_t* E)
{
    ac_reset_decoder(E->ac);
}

