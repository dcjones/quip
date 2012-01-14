
#include "idenc.h"
#include "ac.h"
#include "cumdist.h"
#include "misc.h"
#include <stdlib.h>
#include <string.h>
#include <ctype.h>


struct idenc_t_
{
    ac_t* ac;

    size_t max_id_len;

    /* distribution over edit operations by position */
    cumdist_t** cs;

    /* run lengths of match operations by position */
    cumdist_t** ls;

    /* distribution over new characters by position */
    cumdist_t** ds;

    char* lastid;
    size_t lastid_len, lastid_size;



};


/* edit operations used for delta encoding */
typedef enum {
    EDIT_MAT = 0, /* match   */
    EDIT_REP,     /* replace */
    EDIT_DEL,     /* delete  */
    EDIT_INS      /* insert  */
} edit_t;



idenc_t* idenc_alloc(quip_block_writer_t writer, void* writer_data)
{
    idenc_t* E = malloc_or_die(sizeof(idenc_t));

    E->ac = ac_alloc(writer, writer_data);

    E->max_id_len = 0;

    E->cs = NULL;
    E->ls = NULL;
    E->ds = malloc_or_die(256 * sizeof(cumdist_t*));
    size_t i;
    for (i = 0; i < 256; ++i) {
        E->ds[i] = cumdist_alloc(256);
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

    for (i = 0; i < 4 * E->max_id_len; ++i) {
        cumdist_free(E->cs[i]);
    }
    free(E->cs);

    for (i = 0; i < E->max_id_len; ++i) {
        cumdist_free(E->ls[i]);
    }
    free(E->ls);

    for (i = 0; i < 256; ++i) {
        cumdist_free(E->ds[i]);
    }
    free(E->ds);

    free(E);
}


static void idenc_expand_max_id_len(idenc_t* E)
{
    size_t i;

    E->max_id_len += 1;
    E->cs = realloc_or_die(E->cs, 4 * E->max_id_len * sizeof(cumdist_t*));
    for (i = 0; i < 4; ++i) {
        E->cs[4 * (E->max_id_len - 1) + i] = cumdist_alloc(4);
    }

    E->ls = realloc_or_die(E->ls, E->max_id_len * sizeof(cumdist_t*));
    E->ls[E->max_id_len - 1] = cumdist_alloc(2);
}


void idenc_encode(idenc_t* E, const seq_t* x)
{
    edit_t lastop = 0;
    const str_t* id = &x->id1;
    size_t matlen;
    size_t i, j, k;
    for (i = 0, j = 0; i < id->n + 1;) {
        if (i + 1 > E->max_id_len) idenc_expand_max_id_len(E);

        if (j < E->lastid_len) {
            if (id->s[i] == E->lastid[j]) {

                matlen = 0;
                while (j < E->lastid_len && id->s[i] == E->lastid[j]) {
                    ++matlen;
                    ++i;
                    ++j;

                    if (i + 1 > E->max_id_len) idenc_expand_max_id_len(E);
                }

                ac_update(
                    E->ac, 
                    cumdist_p_norm(E->cs[i * 4 + lastop], EDIT_MAT),
                    cumdist_P_norm(E->cs[i * 4 + lastop], EDIT_MAT));

                if (matlen > cumdist_n(E->ls[i]) - 1) {
                    /* 0s are emitted to encode an expansion of the distribution
                     * */
                    for (k = 0; k < (matlen + 1) - cumdist_n(E->ls[i]); ++k) {
                        ac_update(
                            E->ac, 
                            cumdist_p_norm(E->ls[i], 0),
                            cumdist_P_norm(E->ls[i], 0));
                    }

                    cumdist_expand(E->ls[i], matlen + 1);
                }

                ac_update(
                    E->ac, 
                    cumdist_p_norm(E->ls[i], matlen),
                    cumdist_P_norm(E->ls[i], matlen));

                cumdist_add(E->ls[i], matlen, 1);
                cumdist_add(E->cs[i * 4 + lastop], EDIT_MAT, 1);

                lastop = EDIT_MAT;
            }
            else if (isspace(E->lastid[j])) {
                ac_update(
                    E->ac,
                    cumdist_p_norm(E->cs[i * 4 + lastop], EDIT_INS),
                    cumdist_P_norm(E->cs[i * 4 + lastop], EDIT_INS));

                cumdist_add(E->cs[i * 4 + lastop], EDIT_INS, 1);

                ac_update(
                    E->ac,
                    cumdist_p_norm(E->ds[(int) E->lastid[j]], id->s[i]),
                    cumdist_P_norm(E->ds[(int) E->lastid[j]], id->s[i]));

                cumdist_add(E->ds[(int) E->lastid[j]], id->s[i], 1);

                lastop = EDIT_INS;
                ++i;
            }
            else if(isspace(id->s[i])) {
                ac_update(
                    E->ac,
                    cumdist_p_norm(E->cs[i * 4 + lastop], EDIT_DEL),
                    cumdist_P_norm(E->cs[i * 4 + lastop], EDIT_DEL));

                cumdist_add(E->cs[i * 4 + lastop], EDIT_DEL, 1);

                lastop = EDIT_DEL;
                ++j;
            }
            else {
                ac_update(
                    E->ac,
                    cumdist_p_norm(E->cs[i * 4 + lastop], EDIT_REP),
                    cumdist_P_norm(E->cs[i * 4 + lastop], EDIT_REP));

                cumdist_add(E->cs[i * 4 + lastop], EDIT_REP, 1);

                ac_update(
                    E->ac,
                    cumdist_p_norm(E->ds[(int) E->lastid[j]], id->s[i]),
                    cumdist_P_norm(E->ds[(int) E->lastid[j]], id->s[i]));

                cumdist_add(E->ds[(int) E->lastid[j]], id->s[i], 1);

                lastop = EDIT_REP;
                ++i;
                ++j;
            }
        }
        else {
            ac_update(
                E->ac,
                cumdist_p_norm(E->cs[i * 4 + lastop], EDIT_INS),
                cumdist_P_norm(E->cs[i * 4 + lastop], EDIT_INS));

            cumdist_add(E->cs[i * 4 + lastop], EDIT_INS, 1);

            ac_update(
                E->ac,
                cumdist_p_norm(E->ds[0], id->s[i]),
                cumdist_P_norm(E->ds[0], id->s[i]));

            cumdist_add(E->ds[0], id->s[i], 1);

            lastop = EDIT_INS;
            ++i;
        }
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
    // TODO
}




