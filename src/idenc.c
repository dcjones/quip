
#include "idenc.h"
#include "ac.h"
#include "dist.h"
#include "misc.h"
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>


static const size_t max_group_len  = 100;
static const size_t max_offset     = 100;

/* does a particular character constitute a separator */
static const bool issep[256] =
  { 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

/* groups that do match the previous id are encoded either as a numerical offset
 * or as text. */
typedef enum {
    ID_GROUP_MATCH = 0,
    ID_GROUP_TEXT,
    ID_GROUP_OFF,
    ID_GROUP_NUM
} group_type_t;


struct idenc_t_
{
    ac_t* ac;

    /* distribution over number of leading matches given group */
    dist100_t* ms;

    /* distribution over how the group is stored (numeric offset or text) */
    dist4_t* ts;

    /* distribution over the length text given group */
    dist100_t* ls;

    /* distribution over new characters given group */
    cond_dist128_t* ds;

    /* distribution over numerical offset given group */
    dist100_t* ns;

    /* distribution over the number of bytes used to encode a number */
    dist8_t* bs;

    /* distribution over byte values used when encoding numbers */
    cond_dist256_t* ks;

    /* number of whitespace seperated groups allowed for */
    size_t groups;

    char* lastid;
    size_t lastid_len, lastid_size;

    bool decoder;
};


static void idenc_init(idenc_t* E)
{
    E->groups = 0;

    E->ms = NULL;
    E->ts = NULL;
    E->ls = NULL;
    E->ds = NULL;
    E->ns = NULL;
    E->bs = NULL;
    E->ks = NULL;

    E->lastid = NULL;
    E->lastid_len = 0;
    E->lastid_size = 0;
}



idenc_t* idenc_alloc_encoder(quip_writer_t writer, void* writer_data)
{
    idenc_t* E = malloc_or_die(sizeof(idenc_t));

    E->ac = ac_alloc_encoder(writer, writer_data);

    idenc_init(E);

    return E;
}



idenc_t* idenc_alloc_decoder(quip_reader_t reader, void* reader_data)
{
    idenc_t* E = malloc_or_die(sizeof(idenc_t));

    E->ac = ac_alloc_decoder(reader, reader_data);

    idenc_init(E);

    return E;
}



void idenc_free(idenc_t* E)
{
    free(E->lastid);
    ac_free(E->ac);

    size_t i;
    for (i = 0; i < E->groups; ++i){
        cond_dist128_free(&E->ds[i]);
        cond_dist256_free(&E->ks[i]);
    }

    free(E->ms);
    free(E->ts);
    free(E->ls);
    free(E->ds);
    free(E->ns);
    free(E->bs);
    free(E->ks);

    free(E);
}


static void idenc_add_group(idenc_t* E)
{
    E->groups++;
    E->ms = realloc_or_die(E->ms, E->groups * sizeof(dist100_t));
    E->ts = realloc_or_die(E->ts, E->groups * sizeof(dist4_t));
    E->ls = realloc_or_die(E->ls, E->groups * sizeof(dist100_t));
    E->ds = realloc_or_die(E->ds, E->groups * sizeof(cond_dist128_t));
    E->ns = realloc_or_die(E->ns, E->groups * sizeof(dist100_t));
    E->bs = realloc_or_die(E->bs, E->groups * sizeof(dist8_t));
    E->ks = realloc_or_die(E->ks, E->groups * sizeof(cond_dist256_t));

    dist100_init      (&E->ms[E->groups - 1]);
    dist4_init        (&E->ts[E->groups - 1]);
    dist100_init      (&E->ls[E->groups - 1]);
    cond_dist128_init (&E->ds[E->groups - 1], 128 * 128);
    dist100_init      (&E->ns[E->groups - 1]);
    dist8_init        (&E->bs[E->groups - 1]);
    cond_dist256_init (&E->ks[E->groups - 1], 8);
}



void idenc_encode(idenc_t* E, const seq_t* seq)
{
    const str_t* id = &seq->id1; /* for convenience */

    size_t i; /* position within the current id */
    size_t j; /* position within the last id */

    size_t i0; /* first non-sep position within the current id group */
    size_t j0; /* first non-sep position within the last id group */

    size_t i_last;

    size_t u; /* index of the end of the group in id */
    size_t v; /* index of the end of the group in lost_id */

    size_t matches;    /* number of leading matches */
    size_t group = 0;  /* group number */
    size_t ctx; /* context for 2nd order model of text groups */

    /* are the groups numeric */
    bool is_num, last_is_num;

    /* current and last numeric values */
    unsigned long long x, y;

    for (i = 0, j = 0; i < id->n + 1; ++group) {
        /* make a new group ? */
        if (group >= E->groups) {
            idenc_add_group(E);
        }

        is_num      = true;
        last_is_num = true;

        i_last = i;

        /* consume leading separator matches */
        matches = 0;
        while (i < id->n + 1 &&
               j < E->lastid_len &&
               matches < max_group_len &&
               id->s[i] == E->lastid[j] &&
               issep[(uint8_t) id->s[i]])
        {
            ++i;
            ++j;
            ++matches;
        }
        i0 = i;
        j0 = j;

        /* consume leading non-separator matches */
        while (i < id->n + 1 &&
               j < E->lastid_len &&
               matches < max_group_len &&
               id->s[i] == E->lastid[j] &&
               !issep[(uint8_t) id->s[i]])
        {
            is_num      = is_num && isdigit(id->s[i]);
            last_is_num = last_is_num && isdigit(E->lastid[j]);

            ++i;
            ++j;
            ++matches;
        }

        /* find the end of the current id group */

        /* consume leading non-matching separators when nothing matched */
        u = i;
        is_num = true;
        if (matches == 0) {
            while (u < id->n + 1 &&
                   issep[(uint8_t) id->s[u]] &&
                   (u - i) < max_group_len)
            {
                is_num = false;
                ++u;
            }
        }

        /* consume trailing non-matching characters */
        while (u < id->n + 1 &&
               !issep[(uint8_t) id->s[u]] &&
               matches + (u - i) < max_group_len)
        {
            is_num = is_num && isdigit(id->s[u]);
            ++u;
        }

        /* find the end of the last id group */
        v = j;
        last_is_num = true;
        while (v < E->lastid_len &&
               !issep[(uint8_t) E->lastid[v]] &&
               matches + (v - j) < max_group_len)
        {
            last_is_num = last_is_num && isdigit(E->lastid[v]);
            ++v;
        }

        /* TODO: handle the case when we are at the maximum number of groups */
        /* TODO: handle the case in which we are at the maximum group length */

        /* groups match completely */
        if (i == u && j == v) {
            dist4_encode (E->ac, &E->ts[group], ID_GROUP_MATCH);
            i = u;
            j = v;
            continue;
        }

        /* TODO: handle the case of leading zeros */

        /* difference encoded as a numerical offset */
        if (is_num && last_is_num) {
            x = strtoull(id->s + i0, NULL, 10);
            y = strtoull(E->lastid + j0, NULL, 10);

            if (x > y && x - y < max_offset) {
                dist4_encode   (E->ac, &E->ts[group], ID_GROUP_OFF);
                dist100_encode (E->ac, &E->ms[group], i0 - i_last);
                dist100_encode (E->ac, &E->ns[group], x - y);
                i = u;
                j = v;
                continue;
            }
            else {
                dist4_encode   (E->ac, &E->ts[group], ID_GROUP_NUM);
                dist100_encode (E->ac, &E->ms[group], i0 - i_last);

                size_t bytes = 0;
                y = x;
                while (y > 0) {
                    y >>= 8;
                    ++bytes;
                }

                /* encode the number of bytes needed (i.e., ceil(log2(x))) */
                dist8_encode(E->ac, &E->bs[group], bytes);

                /* encode the bytes */
                while (x > 0) {
                    cond_dist256_encode(E->ac, &E->ks[group], bytes - 1, x & 0xff);
                    x >>= 8;
                    bytes--;
                }

                i = u;
                j = v;
                continue;
            }
        }

        /* difference encoded as text */
        dist4_encode   (E->ac, &E->ts[group], ID_GROUP_TEXT);
        dist100_encode (E->ac, &E->ms[group], matches);
        dist100_encode (E->ac, &E->ls[group], u - i);

        ctx = 0;
        while (i < u) {
            ctx = 128 * (i > 0 ? id->s[i - 1] : 0);
            ctx += j < E->lastid_len ? E->lastid[j] : 0;

            cond_dist128_encode(E->ac, &E->ds[group], ctx, id->s[i]);

            if (j + 1 < v) ++j;
            ++i;
        }

        j = v;
    }


    if (E->lastid_size < id->n + 1) {
        E->lastid_size = id->n + 1;
        free(E->lastid);
        E->lastid = malloc_or_die((id->n + 1) * sizeof(char));
    }
    memcpy(E->lastid, id->s, (id->n + 1) * sizeof(char));
    E->lastid_len = id->n + 1;
}


void idenc_flush(idenc_t* E)
{
    ac_flush_encoder(E->ac);
}


void idenc_decode(idenc_t* E, seq_t* seq)
{
    str_t* id = &seq->id1; /* for convenience */

    group_type_t t;

    size_t i; /* position within the current id */
    size_t j; /* position within the last id */
    size_t j0;

    size_t u;

    size_t matches;    /* number of leading matches */
    size_t group = 0;  /* group number */
    size_t ctx; /* context for 2nd order model of text groups */

    size_t bytes; /* number of bytes used to encode a number */

    /* current and last numeric values */
    unsigned long long x, y, off;

    i = 0;
    j = 0;

    do{
        /* make a new group ? */
        if (group >= E->groups) {
            idenc_add_group(E);
        }

        t = dist4_decode(E->ac, &E->ts[group]);

        if (t == ID_GROUP_MATCH) {
            j0 = j;

            if (E->lastid[j] == '\0') ++j;
            else {
                while (issep[(uint8_t) E->lastid[j]] && j < E->lastid_len) ++j;
                while (!issep[(uint8_t) E->lastid[j]] && j < E->lastid_len) ++j;
            }
            matches = j - j0;

            while (i + matches >= id->size) fastq_expand_str(id);
            memcpy(id->s + i, E->lastid + j0, matches);

            i += matches;

            goto decode_loop_end;
        }

        matches = dist100_decode (E->ac, &E->ms[group]);
        assert(j + matches <= E->lastid_len);

        while (i + matches >= id->size) fastq_expand_str(id);
        memcpy(id->s + i, E->lastid + j, matches);

        i += matches;
        j += matches;

        if (t == ID_GROUP_OFF) {
            off = dist100_decode(E->ac, &E->ns[group]);

            assert(off < max_offset);

            /* this might be a little dangerous */
            y = strtoull(E->lastid + j, NULL, 10);
            x = y + off;

            /* note: 20 = ceil(log10(2**64))  */
            while (i + 20 >= id->size) fastq_expand_str(id);
            i += snprintf(id->s + i, 20, "%llu", x);

            /* consume numbers in lasdid */
            while (j < E->lastid_len && isdigit(E->lastid[j])) ++j;
        }
        else if (t == ID_GROUP_NUM) {
            bytes = dist8_decode(E->ac, &E->bs[group]);

            x = 0;
            size_t b;
            for (b = 0; b < bytes; ++b) {
                x |= cond_dist256_decode(E->ac, &E->ks[group], bytes - b - 1) << (8 * b);
            }

            /* note: 20 = ceil(log10(2**64))  */
            while (i + 20 >= id->size) fastq_expand_str(id);
            i += snprintf(id->s + i, 20, "%llu", x);

            /* consume numbers in lasdid */
            while (j < E->lastid_len && isdigit(E->lastid[j])) ++j;
        }
        else {
            u = i + dist100_decode(E->ac, &E->ls[group]);
            while (u >= id->size) fastq_expand_str(id);

            ctx = 0;
            while (i < u) {
                ctx = 128 * (i > 0 ? id->s[i - 1] : 0);
                ctx += j < E->lastid_len ? E->lastid[j] : 0;

                id->s[i] = (char) cond_dist128_decode(E->ac, &E->ds[group], ctx);

                if (j + 1 < E->lastid_len && !issep[(uint8_t) (j + 1)]) ++j;
                ++i;
            }

            if (j + 1 < E->lastid_len) ++j;
        }

decode_loop_end:

        ++group;

    } while (i == 0 || id->s[i - 1] != '\0');

    id->n = i;


    if (E->lastid_size < id->n + 1) {
        E->lastid_size = id->n + 1;
        free(E->lastid);
        E->lastid = malloc_or_die((id->n + 1) * sizeof(char));
    }
    memcpy(E->lastid, id->s, (id->n + 1) * sizeof(char));
    E->lastid_len = id->n + 1;
}


void idenc_reset_decoder(idenc_t* E)
{
    ac_reset_decoder(E->ac);
}


