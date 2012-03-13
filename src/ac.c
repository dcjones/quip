
#include "ac.h"
#include "misc.h"
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <stdbool.h>

/* allowable values of the interval length before renormalization */
const uint32_t min_length = 0x01000000U;
const uint32_t max_length = 0xFFFFFFFFU;


ac_t* ac_alloc_encoder(quip_writer_t writer, void* writer_data)
{
    ac_t* ac = malloc_or_die(sizeof(ac_t));
    ac->b = 0;
    ac->l = max_length;

    ac->buflen = 10000000;
    ac->buf = malloc_or_die(ac->buflen * sizeof(uint8_t));
    ac->bufpos = 0;

    ac->writer = writer;
    ac->writer_data = writer_data;

    ac->reader = NULL;
    ac->reader_data = NULL;

    return ac;
}


ac_t* ac_alloc_decoder(quip_reader_t reader, void* reader_data)
{
    ac_t* ac = malloc_or_die(sizeof(ac_t));
    ac->l = max_length;

    ac->buflen = 4096;
    ac->buf = malloc_or_die(ac->buflen * sizeof(uint8_t));
    ac->bufpos = 0;

    ac->writer = NULL;
    ac->writer_data = NULL;

    ac->reader = reader;
    ac->reader_data = reader_data;

    return ac;
}


void  ac_free(ac_t* ac)
{
    free(ac->buf);
    free(ac);
}



static void ac_append_byte(ac_t* E, uint8_t c)
{
    if (E->bufpos >= E->buflen) {
        E->buflen += 1000000;
        E->buf = realloc_or_die(E->buf, E->buflen * sizeof(uint8_t));
     }

    assert(E->bufpos < E->buflen);

    E->buf[E->bufpos++] = c;
}


static uint8_t ac_get_byte(ac_t* ac)
{
    if (ac->bufpos < ac->bufavail) {
        return ac->buf[ac->bufpos++];
    }
    else {
        ac->bufavail = ac->reader(ac->reader_data, ac->buf, ac->buflen);
        ac->bufpos = 0;
        if (ac->bufpos < ac->bufavail) return ac->buf[ac->bufpos++];
        else return 0;
    }
}


void ac_propogate_carry(ac_t* ac)
{
    if (ac->bufpos == 0) return;

    size_t i;
    for (i = ac->bufpos - 1; i > 0; --i) {
        if (ac->buf[i] < 0xff) {
            ac->buf[i] += 1;
            break;
        }
        else ac->buf[i] = 0;
    }

    if (i == 0) {
        assert(ac->buf[0] < 0xff);
        ac->buf[0] += 1;
    }
}


void ac_renormalize_encoder(ac_t* ac)
{
    while (ac->l < min_length) {
        ac_append_byte(ac, (uint8_t) (ac->b >> 24));
        ac->b <<= 8;
        ac->l <<= 8;
    }
}


void ac_renormalize_decoder(ac_t* ac)
{
    while (ac->l < min_length) {
        ac->v = (ac->v << 8) | (uint32_t) ac_get_byte(ac);
        ac->l <<= 8;
    }
}


size_t ac_finish_encoder(ac_t* ac)
{
    uint32_t b0 = ac->b;

    if (ac->l > 2 * min_length) {
        ac->b += min_length;
        ac->l = min_length >> 1;
    }
    else {
        ac->b += min_length >> 1;
        ac->l = min_length >> 9;
    }

    if (b0 > ac->b) ac_propogate_carry(ac);
    ac_renormalize_encoder(ac);

    return ac->bufpos;
}


void ac_flush_encoder(ac_t* ac)
{
    ac->writer(ac->writer_data, ac->buf, ac->bufpos);

    ac->b = 0;
    ac->l = max_length;
    ac->bufpos = 0;
}


void ac_start_decoder(ac_t* ac)
{
    ac->bufavail = ac->reader(ac->reader_data, ac->buf, ac->buflen);

    if (ac->bufavail < 4) {
        fprintf(stderr, "Malformed compressed data encountered.");
        exit(EXIT_FAILURE);
    }

    ac->v = ((uint32_t) ac->buf[0] << 24) | ((uint32_t) ac->buf[1] << 16) |
    ((uint32_t) ac->buf[2] << 8)  | ((uint32_t) ac->buf[3]);

    ac->l = max_length;

    ac->bufpos = 4;
}


void ac_reset_decoder(ac_t* ac)
{
    ac->bufpos = 0;
    ac->bufavail = 0;
}


