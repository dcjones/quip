
#include "idenc.h"
#include "ac.h"
#include "dist.h"
#include "misc.h"
#include <stdlib.h>
#include <string.h>
#include <ctype.h>


/* How many dist_add calls are made before all the distributions are updated. */
static const size_t dist_update_delay = 10000;

static const size_t max_group_len = 100;
static const size_t max_groups = 10;

struct idenc_t_
{
    ac_t* ac;

    size_t max_id_len;

    /* distribution over number of leading matches by whitespace seperated group
     * */
    dist_t** ms;

    /* distribution over new characters by group */
    dist_t** ds;

    /* number of whitespace seperated groups allowed for */
    size_t groups;

    char* lastid;
    size_t lastid_len, lastid_size;
};


static void idenc_init(idenc_t* E, bool decode)
{
    E->max_id_len = 0;

    E->groups = 1;

    size_t i;

    E->ms = malloc_or_die(E->groups * sizeof(dist_t*));
    for (i = 0; i < E->groups; ++i) {
        E->ms[i] = dist_alloc(max_group_len, decode);
    }

    /* TODO: allocate these lazily to conserve space */
    E->ds = malloc_or_die(max_groups * 128 * 128 * sizeof(dist_t*));
    for (i = 0; i < max_groups * 128 * 128; ++i) {
        E->ds[i] = dist_alloc(128, decode);
    }

    E->lastid = NULL;
    E->lastid_len = 0;
    E->lastid_size = 0;
}



idenc_t* idenc_alloc_encoder(quip_writer_t writer, void* writer_data)
{
    idenc_t* E = malloc_or_die(sizeof(idenc_t));

    E->ac = ac_alloc_encoder(writer, writer_data);

    idenc_init(E, false);

    return E;
}



idenc_t* idenc_alloc_decoder(quip_reader_t reader, void* reader_data)
{
    idenc_t* E = malloc_or_die(sizeof(idenc_t));

    E->ac = ac_alloc_decoder(reader, reader_data);

    idenc_init(E, true);

    return E;
}



void idenc_free(idenc_t* E)
{
    free(E->lastid);
    ac_free(E->ac);
    size_t i;
    for (i = 0; i < E->groups; ++i) {
        dist_free(E->ms[i]);
    }
    free(E->ms);

    for (i = 0; i < max_groups * 128 * 128; ++i) {
        dist_free(E->ds[i]);
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
    size_t i, j;
    size_t matches;
    size_t group = 0;

    size_t ctx;

    const str_t* id = &x->id1;

    for (i = 0, j = 0; i < id->n + 1; ++group) {
        /* make a new group ? */
        if (group >= E->groups) {
            E->groups++;
            E->ms = realloc_or_die(E->ms, E->groups * sizeof(dist_t*));
            E->ms[group] = dist_alloc_encode(max_group_len);
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

        ac_encode(E->ac, E->ms[group], matches);

        /* write trailing whitespace */
        while (i < id->n + 1 && issep(id->s[i])) {
            ctx = 128 * 128 * (group >= max_groups ? max_groups - 1 : group);
            ctx += 128 * (i > 0 ? id->s[i - 1] : 0);
            /*c += 128 * (j - 1 < E->lastid_len ? E->lastid[j - 1] : 0);*/
            ctx += j < E->lastid_len ? E->lastid[j] : 0;
            ctx %= max_groups * 128 * 128;

            ac_encode(E->ac, E->ds[ctx], id->s[i]);

            ++i;

            if (j < E->lastid_len && issep(E->lastid[j])) ++j;
        }

        /* write non-matches */
        while (i < id->n + 1 && !issep(id->s[i])) {
            ctx = 128 * 128 * (group >= max_groups ? max_groups - 1 : group);
            ctx += 128 * (i > 0 ? id->s[i - 1] : 0);
            /*c += 128 * (j - 1 < E->lastid_len ? E->lastid[j - 1] : 0);*/
            ctx += j < E->lastid_len ? E->lastid[j] : 0;
            ctx %= max_groups * 128 * 128;

            ac_encode(E->ac, E->ds[ctx], id->s[i]);
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
    ac_flush_encoder(E->ac);
}

// TODO:
/* The below is a fucking distaster. */
void idenc_decode(idenc_t* E, seq_t* seq)
{
    size_t i = 0;
    size_t j = 0;
    size_t group = 0;
    size_t matches;

    size_t ctx;
    uint8_t x;

    str_t* id = &seq->id1;

    bool consume_sep;
    
    while (true) {
        /* make a new group ? */
        if (group >= E->groups) {
            E->groups++;
            E->ms = realloc_or_die(E->ms, E->groups * sizeof(dist_t*));
            E->ms[group] = dist_alloc_encode(max_group_len);
        }

        matches = (size_t) ac_decode(E->ac, E->ms[group]);
        while (i + matches >= id->size) fastq_expand_str(id);
        memcpy(id->s + i, E->lastid + j, matches);
        i += matches;
        j += matches;

        do {
            ctx = 128 * 128 * (group >= max_groups ? max_groups - 1 : group);
            ctx += 128 * (i > 0 ? id->s[i - 1] : 0);
            ctx += j < E->lastid_len ? E->lastid[j] : 0;
            ctx %= max_groups * 128 * 128;

            x = ac_decode(E->ac, E->ds[ctx]);
            while (id->n + 1 >= id->size) fastq_expand_str(id);
            id->s[i++] = x;

            if (j < E->lastid_len && !issep(E->lastid[j])) ++j;


        } while (true); // TODO

        /* scan through lastid group */
        while (j < E->lastid_len && !issep(E->lastid[j])) ++j;
    }

    id->s[i] = '\0';
    id->n = i;
}




