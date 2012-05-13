
#include "fastqfmt.h"
#include "misc.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

struct quip_fastq_out_t_
{
    quip_writer_t writer;
    void*         writer_data;
    char*         buf;
    size_t        buflen;
};


quip_fastq_out_t* quip_fastq_out_open(
                    quip_writer_t writer,
                    void*         writer_data,
                    ATTRIB_UNUSED quip_opt_t opts)
{
    quip_fastq_out_t* out = malloc_or_die(sizeof(quip_fastq_out_t));
    out->writer      = writer;
    out->writer_data = writer_data;
    out->buflen      = 1024;
    out->buf         = malloc_or_die(out->buflen);
    return out;
}


void quip_fastq_out_close(quip_fastq_out_t* out)
{
    if (out) {
        free(out->buf);
        free(out);
    }
}


void quip_fastq_write(quip_fastq_out_t* out, short_read_t* sr)
{
    size_t size_needed = sr->id.n + sr->seq.n + sr->qual.n + 7;

    if (size_needed > out->buflen) {
        out->buflen = size_needed;
        out->buf = realloc_or_die(out->buf, out->buflen * sizeof(char));
    }

    snprintf(out->buf, size_needed, "@%s\n%s\n+\n%s\n", sr->id.s, sr->seq.s, sr->qual.s);
    out->writer(out->writer_data, (uint8_t*) out->buf, size_needed - 1);
}



static const size_t fastq_in_bufsize = 1048576;

typedef enum
{
    STATE_EOF,
    STATE_ERROR,
    STATE_ID1,
    STATE_SEQ,
    STATE_ID2,
    STATE_QUAL

} fastq_state_t;


struct quip_fastq_in_t_
{
    quip_reader_t reader;
    void*         reader_data;
    fastq_state_t state;
    char* buf;
    char* c;
    short_read_t sr;
};

quip_fastq_in_t* quip_fastq_in_open(
                    quip_reader_t reader,
                    void*         reader_data,
                    ATTRIB_UNUSED quip_opt_t opts)
{
    quip_fastq_in_t* in = malloc_or_die(sizeof(quip_fastq_in_t));

    in->reader      = reader;
    in->reader_data = reader_data;
    in->state       = STATE_ID1;
    in->buf         = malloc_or_die(fastq_in_bufsize);
    in->buf[0]      = '\0';
    in->c           = in->buf;
    short_read_init(&in->sr);

    return in;
}


void quip_fastq_in_close(quip_fastq_in_t* in)
{
    if (in) {
        free(in->buf);
        short_read_free(&in->sr);
        free(in);
    }
}


static void fastq_in_refill(quip_fastq_in_t* in)
{
    size_t n = in->reader(
        in->reader_data, (uint8_t*) in->buf, fastq_in_bufsize - 1);

    if (n == 0) {
        in->state = STATE_EOF;
    }
    else {
        in->buf[n] = '\0';
        in->c = in->buf;
    }
}


static void fastq_in_get_line(quip_fastq_in_t* in, str_t* s)
{
    if (s) s->n = 0;

    size_t i = 0;
    while (1) {
        while (in->c[i] >= '\r') ++i;

        if (in->c[i] == '\0') {
            if (s)  {
                str_reserve(s, s->n + i + 2);
                memcpy(s->s + s->n, in->c, i);
                s->n += i;
            }

            fastq_in_refill(in);
            i = 0;
            if (in->state == STATE_EOF) break;
        }
        else {
            if (s)  {
                str_reserve(s, s->n + i + 2);
                memcpy(s->s + s->n, in->c, i);
                s->n += i;
            }

            in->c += i;

            break;
        }
    }

    if (s) {
        if (s->n > 0 && s->s[s->n - 1] == '\r') s->n--;
        s->s[s->n] = '\0';
    }
}


short_read_t* quip_fastq_read(quip_fastq_in_t* in)
{
    if (in->state == STATE_EOF) return NULL;

    while (true) {
        /* read more, if needed */
        if (*in->c == '\0' ) {
            fastq_in_refill(in);
            if (in->state == STATE_EOF) return NULL;
            continue;
        }
        else if (isspace(*in->c)) {
            in->c++;
            continue;
        }

        switch (in->state) {
            case STATE_ID1:
                if (*in->c == '@' || *in->c == '>') {
                    in->c++;
                    fastq_in_get_line(in, &in->sr.id);
                    if (in->state == STATE_EOF) return NULL;

                    in->state = STATE_SEQ;
                }
                else {
                    quip_error(
                        "Malformed FASTQ file: expecting an '@' or '>', saw a '%c'");
                }
                break;


            case STATE_SEQ:
                fastq_in_get_line(in, &in->sr.seq);
                if (in->state == STATE_EOF) return NULL;

                in->state = STATE_ID2;
                break;


            case STATE_ID2:
                if (*in->c == '+') {
                    in->c++;
                    fastq_in_get_line(in, NULL);
                    if (in->state == STATE_EOF) return NULL;

                    in->state = STATE_QUAL;
                }
                else {
                    /* fasta style entry */
                    if (in->sr.qual.size > 0) {
                        in->sr.qual.s[0] = '\0';
                        in->sr.qual.n = 0;
                    }

                    in->state = STATE_ID1;
                    return &in->sr;
                }
                break;


            case STATE_QUAL:
                fastq_in_get_line(in, &in->sr.qual);
                if (in->state == STATE_EOF) return &in->sr;

                in->state = STATE_ID1;
                return &in->sr;

            default: break;
        }

        in->c++;
    }

    return &in->sr;
}

