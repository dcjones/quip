
#include "idenc.h"
#include "ac.h"
#include "dist.h"
#include "misc.h"
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>


static const size_t max_group_len  = 100;
static const size_t max_groups     = 20;
static const size_t max_offset     = 100;

/* groups that do match the previous id are encoded either as a numerical offset
 * or as text. */
typedef enum {
    ID_GROUP_MATCH = 0,
    ID_GROUP_TEXT,
    ID_GROUP_NUM
} group_type_t;


struct idenc_t_
{
    ac_t* ac;

    size_t max_id_len;

    /* distribution over number of leading matches given group */
    dist_t** ms;

    /* distribution over how the group is stored (numeric offset or text) */
    dist_t** ts;

    /* distribution over the length text given group */
    dist_t** ls;

    /* distribution over new characters given group */
    dist_t*** ds;

    /* distribution over numerical offset given group */
    dist_t** ns;

    /* number of whitespace seperated groups allowed for */
    size_t groups;

    char* lastid;
    size_t lastid_len, lastid_size;

    bool decoder;
};


static void idenc_init(idenc_t* E)
{
    E->max_id_len = 0;
    E->groups = 0;

    E->ms = malloc_or_die(max_groups * sizeof(dist_t*));
    memset(E->ms, 0, max_groups * sizeof(dist_t*));

    E->ts = malloc_or_die(max_groups * sizeof(dist_t*));
    memset(E->ts, 0, max_groups * sizeof(dist_t*));

    E->ls = malloc_or_die(max_groups * sizeof(dist_t*));
    memset(E->ls, 0, max_groups * sizeof(dist_t*));

    E->ds = malloc_or_die(max_groups * sizeof(dist_t**));
    memset(E->ds, 0, max_groups * sizeof(dist_t**));

    E->ns = malloc_or_die(max_groups * sizeof(dist_t*));
    memset(E->ns, 0, max_groups * sizeof(dist_t*));

    E->lastid = NULL;
    E->lastid_len = 0;
    E->lastid_size = 0;
}



idenc_t* idenc_alloc_encoder(quip_writer_t writer, void* writer_data)
{
    idenc_t* E = malloc_or_die(sizeof(idenc_t));

    E->ac = ac_alloc_encoder(writer, writer_data);

    idenc_init(E);

    E->decoder = false;

    return E;
}



idenc_t* idenc_alloc_decoder(quip_reader_t reader, void* reader_data)
{
    idenc_t* E = malloc_or_die(sizeof(idenc_t));

    E->ac = ac_alloc_decoder(reader, reader_data);

    idenc_init(E);

    E->decoder = true;

    return E;
}



void idenc_free(idenc_t* E)
{
    free(E->lastid);
    ac_free(E->ac);

    size_t i, j;
    for (i = 0; i < max_groups; ++i) dist_free(E->ms[i]);
    free(E->ms);

    for (i = 0; i < max_groups; ++i) dist_free(E->ts[i]);
    free(E->ts);

    for (i = 0; i < max_groups; ++i) dist_free(E->ls[i]);
    free(E->ls);

    for (i = 0; i < max_groups; ++i) {
        if (E->ds[i]) {
            for (j = 0; j < 128 * 128; ++j) {
                dist_free(E->ds[i][j]);
            }
        }
        free(E->ds[i]);
    }
    free(E->ds);

    for (i = 0; i < max_groups; ++i) dist_free(E->ns[i]);
    free(E->ns);

    free(E);
}


static bool issep(char c)
{
    return isspace(c) || c == '/' || c == '.' || c == ':';
}


static bool idenc_add_group(idenc_t* E)
{
    if (E->groups == max_groups) return false;

    E->ms[E->groups] = dist_alloc(max_group_len, E->decoder);
    E->ts[E->groups] = dist_alloc(3, E->decoder);
    E->ls[E->groups] = dist_alloc(max_group_len, E->decoder);

    E->ds[E->groups] = malloc_or_die(128 * 128 * sizeof(dist_t*));
    size_t j;
    for (j = 0; j < 128 * 128; ++j) {
        E->ds[E->groups][j] = dist_alloc(128, E->decoder);
    }

    E->ns[E->groups] = dist_alloc(max_offset + 1, E->decoder);

    E->groups++;

    return true;
}



void idenc_encode(idenc_t* E, const seq_t* seq)
{
    const str_t* id = &seq->id1; /* for convenience */

    size_t i; /* position within the current id */
    size_t j; /* position within the last id */

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
            assert(group < max_groups - 1);
            idenc_add_group(E);
        }

        is_num      = true;
        last_is_num = true;

        /* consume leading separator matches */
        matches = 0;
        while (i < id->n + 1 &&
               j < E->lastid_len &&
               matches < max_group_len &&
               id->s[i] == E->lastid[j] &&
               issep(id->s[i]))
        {
            is_num      = false;
            last_is_num = false;
            is_num      = is_num && isdigit(id->s[i]);
            last_is_num = last_is_num && isdigit(E->lastid[j]);

            ++i;
            ++j;
            ++matches;
        }


        /* consume leading non-separator matches */
        while (i < id->n + 1 &&
               j < E->lastid_len &&
               matches < max_group_len &&
               id->s[i] == E->lastid[j] &&
               !issep(id->s[i]))
        {
            is_num      = is_num && isdigit(id->s[i]);
            last_is_num = last_is_num && isdigit(E->lastid[j]);

            ++i;
            ++j;
            ++matches;
        }


        /* find trailing numbers */
        if (is_num && last_is_num) {
            u = i;
            while (u <- id->n + 1 &&
                   isdigit(id->s[u]) &&
                   matches + (u - i) < max_group_len)
            {
                ++u;
            }

            v = j;
            while (v <- E->lastid_len &&
                   isdigit(E->lastid[v]) &&
                   matches + (v - j) < max_group_len)
            {
                ++v;
            }

            if (matches > 0 || ((u - i) > 0 && (v - j) > 0)) {
                x = strtoull(id->s + (i - matches), NULL, 10);
                y = strtoull(E->lastid + (j - matches), NULL, 10);

                if (x > y && x - y <= max_offset) {
                    ac_encode(E->ac, E->ts[group], ID_GROUP_NUM);
                    ac_encode(E->ac, E->ns[group], x - y);

                    i = u;
                    j = v;

                    continue;
                }
            }
        }


        /* encode matches */
        if (matches > 0) {
            /* try to consume the trailing seperator */
            /*if (i < id->n + 1 &&*/
                /*j < E->lastid_len &&*/
                /*matches < max_group_len &&*/
                /*id->s[i] == E->lastid[j]) {*/
                /*++i;*/
                /*++j;*/
                /*++matches;*/
            /*}*/

            ac_encode(E->ac, E->ts[group], ID_GROUP_MATCH);
            ac_encode(E->ac, E->ms[group], matches);
            continue;
        }


        /* encode non-matches up to and including the first separator */

        /* find the end of the current id group */

        /* consume leading non-matching separators when nothing matched */
        u = i;
        while (u < id->n + 1 &&
               issep(id->s[u]) &&
               matches + (u - i) < max_group_len)
        {
            ++u;
        }

        /* consume trailing non-matching characters */
        while (u < id->n + 1 &&
               !issep(id->s[u]) &&
               matches + (u - i) < max_group_len)
        {
            ++u;
        }

        /* include the final separator */
        /*if (u < id->n + 1 && matches + (u - i) < max_group_len) ++u;*/

        /* find the end of the last id group */
        v = j;
        last_is_num = true;
        while (v < E->lastid_len &&
               !issep(E->lastid[v]) &&
               matches + (v - j) < max_group_len)
        {
            ++v;
        }

        /* include the final separator */
        /*if (v < E->lastid_len && matches + (v - j) < max_group_len) ++v;*/

        /* TODO: handle the case when we are at the maximum number of groups */
        /* TODO: handle the case in which we are at the maximum group length */

        ac_encode(E->ac, E->ts[group], ID_GROUP_TEXT);
        ac_encode(E->ac, E->ls[group], u - i);

        ctx = 0;
        while (i < u) {
            ctx = 128 * (i > 0 ? id->s[i - 1] : 0);
            ctx += j < E->lastid_len ? E->lastid[j] : 0;

            ac_encode(E->ac, E->ds[group][ctx], id->s[i]);

            if (j + 1 < v) ++j;
            ++i;
        }

        j = v;



#if 0
        /* groups match completely */
        if (i == u && j == v) {
            ac_encode(E->ac, E->ts[group], ID_GROUP_MATCH);
            ac_encode(E->ac, E->ms[group], matches);
            i = u;
            j = v;
            continue;
        }

        /* difference encoded as a numerical offset */
        if (is_num && last_is_num) {
            x = strtoull(id->s + i0, NULL, 10);
            y = strtoull(E->lastid + j0, NULL, 10);

            if (x > y && x - y <= max_offset) {
                ac_encode(E->ac, E->ts[group], ID_GROUP_NUM);
                ac_encode(E->ac, E->ms[group], matches);
                ac_encode(E->ac, E->ns[group], x - y);
                i = u;
                j = v;
                continue;
            }
        }

        /* difference encoded as text */
        ac_encode(E->ac, E->ts[group], ID_GROUP_TEXT);
        ac_encode(E->ac, E->ms[group], matches);
        ac_encode(E->ac, E->ls[group], u - i);

        ctx = 0;
        while (i < u) {
            ctx = 128 * (i > 0 ? id->s[i - 1] : 0);
            ctx += j < E->lastid_len ? E->lastid[j] : 0;

            ac_encode(E->ac, E->ds[group][ctx], id->s[i]);

            if (j + 1 < v) ++j;
            ++i;
        }

        j = v;
#endif
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
    /* TODO */
}




