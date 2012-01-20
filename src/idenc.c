
#include "idenc.h"
#include "ac.h"
#include "cumdist.h"
#include "misc.h"
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

const size_t max_group_len = 100;
const size_t max_groups = 10;

struct idenc_t_
{
    ac_t* ac;

    size_t max_id_len;

    /* distribution over number of leading matches by whitespace seperated group
     * */
    cumdist_t** ms;

    /* distribution over new characters by group */
    cumdist_t** ds;

    /* number of whitespace seperated groups allowed for */
    size_t groups;

    char* lastid;
    size_t lastid_len, lastid_size;
};


idenc_t* idenc_alloc(quip_block_writer_t writer, void* writer_data)
{
    idenc_t* E = malloc_or_die(sizeof(idenc_t));

    E->ac = ac_alloc(writer, writer_data);

    E->max_id_len = 0;

    E->groups = 1;

    size_t i;

    E->ms = malloc_or_die(E->groups * sizeof(cumdist_t*));
    for (i = 0; i < E->groups; ++i) {
        E->ms[i] = cumdist_alloc(max_group_len);
    }

    E->ds = malloc_or_die(max_groups * 128 * 128 * sizeof(cumdist_t*));
    for (i = 0; i < max_groups * 128 * 128; ++i) {
        E->ds[i] = cumdist_alloc(128);
    }

    E->lastid = NULL;
    E->lastid_len = 0;
    E->lastid_size = 0;

    return E;
}



void idenc_free(idenc_t* E)
{
    free(E->lastid);
    ac_free(E->ac);
    size_t i;
    for (i = 0; i < E->groups; ++i) {
        cumdist_free(E->ms[i]);
    }
    free(E->ms);

    for (i = 0; i < max_groups * 128 * 128; ++i) {
        cumdist_free(E->ds[i]);
    }

    free(E->ds);

    free(E);
}


static bool issep(char c)
{
    return isspace(c) || c == '/';
}


void idenc_encode(idenc_t* E, const seq_t* x)
{
    uint32_t p, c;
    size_t i, j;
    size_t matches;
    size_t group = 0;

    size_t ctx;

    const str_t* id = &x->id1;

    for (i = 0, j = 0; i < id->n + 1; ++group) {
        /* make a new group ? */
        if (group >= E->groups) {
            E->groups++;
            E->ms = realloc_or_die(E->ms, E->groups * sizeof(cumdist_t*));
            E->ms[group] = cumdist_alloc(max_group_len);
        }

        /* consume matching whitespace */
        matches = 0;
        while (i < id->n + 1 && j < E->lastid_len &&
               issep(id->s[i]) && issep(E->lastid[j]) &&
               matches < max_group_len && id->s[i] == E->lastid[j])
        {
            ++i;
            ++j;
            ++matches;
        }

        /* consume matching non-whitespace */
        while (i < id->n + 1 && j < E->lastid_len &&
               !issep(id->s[i]) && !issep(E->lastid[j]) &&
               matches < max_group_len && id->s[i] == E->lastid[j])
        {
            ++i;
            ++j;
            ++matches;
        }

        p = cumdist_p_norm(E->ms[group], matches);
        c = cumdist_c_norm(E->ms[group], matches);
        ac_update(E->ac, p, c);
        cumdist_add(E->ms[group], matches, 1);

        /* write trailing whitespace */
        while (i < id->n + 1 && issep(id->s[i])) {
            ctx = 128 * 128 * (group >= max_groups ? max_groups - 1 : group);
            ctx += 128 * (i - 1 < id->n ? id->s[i - 1] : 0);
            /*c += 128 * (j - 1 < E->lastid_len ? E->lastid[j - 1] : 0);*/
            ctx += j < E->lastid_len ? E->lastid[j] : 0;
            ctx %= max_groups * 128 * 128;

            p = cumdist_p_norm(E->ds[ctx], id->s[i]);
            c = cumdist_c_norm(E->ds[ctx], id->s[i]);
            ac_update(E->ac, p, c);
            cumdist_add(E->ds[ctx], id->s[i], 1);
            ++i;

            if (j < E->lastid_len && issep(E->lastid[j])) ++j;
        }

        /* write non-matches */
        while (i < id->n + 1 && !issep(id->s[i])) {
            ctx = 128 * 128 * (group >= max_groups ? max_groups - 1 : group);
            ctx += 128 * (i - 1 < id->n ? id->s[i - 1] : 0);
            /*c += 128 * (j - 1 < E->lastid_len ? E->lastid[j - 1] : 0);*/
            ctx += j < E->lastid_len ? E->lastid[j] : 0;
            ctx %= max_groups * 128 * 128;

            p = cumdist_p_norm(E->ds[ctx], id->s[i]);
            c = cumdist_c_norm(E->ds[ctx], id->s[i]);
            ac_update(E->ac, p, c);
            cumdist_add(E->ds[ctx], id->s[i], 1);
            ++i;

            if (j < E->lastid_len && !issep(E->lastid[j])) ++j;
        }


        /* scan through lastid group */
        while (j < E->lastid_len && !issep(E->lastid[j])) ++j;
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
    ac_flush(E->ac);
}




