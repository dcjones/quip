
#include "idenc.h"
#include "misc.h"
#include "lzma/api/lzma.h"
#include <stdlib.h>
#include <string.h>


struct idenc_t_
{
    quip_block_writer_t writer;
    void* writer_data;

    lzma_stream strm;

    uint8_t out_buf[4096];
};


idenc_t* idenc_alloc(quip_block_writer_t writer, void* writer_data)
{
    idenc_t* E = malloc_or_die(sizeof(idenc_t));

    E->writer = writer;
    E->writer_data = writer_data;

    memset(&E->strm, 0, sizeof(lzma_stream));

    lzma_ret ret = lzma_easy_encoder(&E->strm, 6 /* compression level in [0, 9] */, LZMA_CHECK_CRC32);
    if (ret != LZMA_OK) {
        fprintf(stderr, "lzma_easy_encoder error: %d\n", (int) ret);
        exit(EXIT_FAILURE);
    }

    return E;
}

void idenc_free(idenc_t* E)
{
    lzma_end(&E->strm);
    free(E);
}

void idenc_encode(idenc_t* E, const seq_t* x)
{
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
}


void idenc_flush(idenc_t* E)
{
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
}




