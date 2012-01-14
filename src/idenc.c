
#include "idenc.h"
#include "ac.h"
#include "cumdist.h"
#include "misc.h"
#include "lzma/api/lzma.h"
#include <stdlib.h>
#include <string.h>
#include <ctype.h>


struct idenc_t_
{
    quip_block_writer_t writer;
    void* writer_data;

    lzma_stream strm;

    uint8_t out_buf[4096];


    // NEW MODEL
    ac_t* ac;

    size_t max_id_len;

    /* distribution over edit operations by position */
    cumdist_t** cs;

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

    E->writer = writer;
    E->writer_data = writer_data;

    memset(&E->strm, 0, sizeof(lzma_stream));

    lzma_ret ret = lzma_easy_encoder(&E->strm, 9 /* compression level in [0, 9] */, LZMA_CHECK_CRC32);
    if (ret != LZMA_OK) {
        fprintf(stderr, "lzma_easy_encoder error: %d\n", (int) ret);
        exit(EXIT_FAILURE);
    }

    E->ac = ac_alloc(writer, writer_data);

    E->max_id_len = 0;

    E->cs = NULL;
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
    lzma_end(&E->strm);
    free(E->lastid);
    ac_free(E->ac);
    size_t i;
    for (i = 0; i < 4 * E->max_id_len; ++i) {
        cumdist_free(E->cs[i]);
    }
    free(E->cs);
    for (i = 0; i < 256; ++i) {
        cumdist_free(E->ds[i]);
    }
    free(E->ds);
    free(E);
}

void idenc_encode(idenc_t* E, const seq_t* x)
{
    edit_t lastop = 0;
    const str_t* id = &x->id1;
    size_t i, j, k;
    for (i = 0, j = 0; i < id->n + 1;) {
        if (i + 1 > E->max_id_len) {
            E->max_id_len += 1;
            E->cs = realloc_or_die(E->cs, 4 * E->max_id_len * sizeof(cumdist_t*));
            for (k = 0; k < 4; ++k) {
                E->cs[4 * (E->max_id_len - 1) + k] = cumdist_alloc(4);
            }
        }

        if (j < E->lastid_len) {
            if (id->s[i] == E->lastid[j]) {
                ac_update(
                    E->ac, 
                    cumdist_p_norm(E->cs[i * 4 + lastop], EDIT_MAT),
                    cumdist_P_norm(E->cs[i * 4 + lastop], EDIT_MAT));

                cumdist_add(E->cs[i * 4 + lastop], EDIT_MAT, 1);

                lastop = EDIT_MAT;
                ++i;
                ++j;
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

#if 0
    E->strm.next_in  = (uint8_t*) x->id1.s;
    E->strm.avail_in = x->id1.n + 1; /* +1 to include the trailing '\0' byte */

    lzma_ret ret;

    do {
        E->strm.next_out  = E->out_buf;
        E->strm.avail_out = sizeof(E->out_buf);

        ret = lzma_code(&E->strm, LZMA_RUN);

        if (ret != LZMA_OK && ret != LZMA_STREAM_END) {
            fprintf(stderr, "lzma_code error: %d\n", (int) ret);
            exit(EXIT_FAILURE);
        }
        else if (E->strm.avail_out < sizeof(E->out_buf)) {
            E->writer(E->writer_data,
                      E->out_buf,
                      sizeof(E->out_buf) - E->strm.avail_out);
        }
    } while (E->strm.avail_out == 0);
#endif
}


void idenc_flush(idenc_t* E)
{
#if 0
    lzma_ret ret;

    do {
        E->strm.next_out  = E->out_buf;
        E->strm.avail_out = sizeof(E->out_buf);

        ret = lzma_code(&E->strm, LZMA_FULL_FLUSH);

        if (ret != LZMA_OK && ret != LZMA_STREAM_END) {
            fprintf(stderr, "lzma_code error: %d\n", (int) ret);
            exit(EXIT_FAILURE);
        }
        else if (E->strm.avail_out < sizeof(E->out_buf)) {
            E->writer(E->writer_data,
                      E->out_buf,
                      sizeof(E->out_buf) - E->strm.avail_out);
        }

    } while (ret != LZMA_STREAM_END);
#endif
}




