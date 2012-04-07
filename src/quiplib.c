
#include "quip.h"
#include "misc.h"
#include "quipfmt.h"
#include "fastqfmt.h"
#include "samfmt.h"
#include "sam/bam.h"
#include <stdlib.h>

bool quip_verbose = false;


void str_init(str_t* str)
{
    str->n    = 0;
    str->size = 0;
    str->s    = NULL;
}

void str_reserve(str_t* str, size_t size)
{
    if (size > str->size) {
        str->size = size;
        str->s = realloc_or_die(str->s, str->size * sizeof(char));
    }
}

void str_free(str_t* str)
{
    if (str) free(str->s);
}


void str_copy(str_t* dest, const str_t* src)
{
    str_reserve(dest, src->n + 1);
    memcpy(dest->s, src->s, (src->n + 1) * sizeof(unsigned char));
    dest->n = src->n;
}


void cigar_init(cigar_t* cig)
{
    cig->n    = 0;
    cig->size = 0;
    cig->ops  = NULL;
    cig->lens = NULL;
}


void cigar_reserve(cigar_t* cig, size_t size)
{
    if (size > cig->size) {
        cig->size = size;
        cig->ops  = realloc_or_die(cig->ops,  cig->size * sizeof(uint8_t));
        cig->lens = realloc_or_die(cig->lens, cig->size * sizeof(uint32_t));
    }
}


void cigar_free(cigar_t* cig)
{
    if (cig) {
        free(cig->ops);
        free(cig->lens);
    }
}


void short_read_init(short_read_t* sr)
{
    str_init(&sr->id);
    str_init(&sr->seq);
    str_init(&sr->qual);
    sr->flags = BAM_FUNMAP;
    str_init(&sr->seqname);
    sr->strand = 0;
    sr->pos    = 0;
    str_init(&sr->aux);
}


void short_read_free(short_read_t* sr)
{
    if (sr) {
        str_free(&sr->id);
        str_free(&sr->seq);
        str_free(&sr->qual);
        str_free(&sr->seqname);
        str_free(&sr->aux);
        free(sr);
    }
}


void short_read_copy(short_read_t* dest, const short_read_t* src)
{
    str_copy(&dest->id,   &src->id);
    str_copy(&dest->seq,  &src->seq);
    str_copy(&dest->qual, &src->qual);
    dest->flags = src->flags;
    str_copy(&dest->seqname, &src->seqname);
    dest->strand = src->strand;
    dest->pos    = src->pos;
    str_copy(&dest->aux, &src->aux);
}


struct quip_out_t_
{
    quip_fmt_t fmt;
    union {
        quip_fastq_out_t* fastq;
        quip_sam_out_t*   sam;
        quip_quip_out_t*  quip;
    } x;
};


quip_out_t* quip_out_open(
              quip_writer_t writer,
              void*         writer_data,
              quip_fmt_t    fmt,
              quip_opt_t    opts)
{
    quip_out_t* out = malloc_or_die(sizeof(quip_out_t));
    out->fmt = fmt;

    switch (fmt) {
        case QUIP_FMT_FASTQ:
            out->x.fastq = quip_fastq_out_open(writer, writer_data, opts);
            break;

        case QUIP_FMT_BAM:
            opts |= QUIP_OPT_SAM_BAM;

        case QUIP_FMT_SAM:
            out->x.sam = quip_sam_out_open(writer, writer_data, opts);
            break;

        case QUIP_FMT_QUIP:
            out->x.quip = quip_quip_out_open(writer, writer_data, opts);
            break;

        case QUIP_FMT_UNDEFINED:
            fprintf(stderr, "Undefined format given.\n");
            exit(EXIT_FAILURE);

        case QUIP_FMT_NULL:
            fprintf(stderr, "Null format given.\n");
            exit(EXIT_FAILURE);
    }

    return out;
}

void quip_out_close(quip_out_t* out)
{
    if (out == NULL) return;

    switch (out->fmt) {
        case QUIP_FMT_FASTQ:
            if (out->x.fastq) quip_fastq_out_close(out->x.fastq);
            break;

        case QUIP_FMT_BAM:
        case QUIP_FMT_SAM:
            if (out->x.sam) quip_sam_out_close(out->x.sam);
            break;

        case QUIP_FMT_QUIP:
            if (out->x.quip) quip_quip_out_close(out->x.quip);
            break;

        default: break;
    }

    free(out);
}


void quip_write(quip_out_t* out, short_read_t* sr)
{
    switch (out->fmt) {
        case QUIP_FMT_FASTQ:
            quip_fastq_write(out->x.fastq, sr);
            break;

        case QUIP_FMT_BAM:
        case QUIP_FMT_SAM:
            quip_sam_write(out->x.sam, sr);
            break;

        case QUIP_FMT_QUIP:
            quip_quip_write(out->x.quip, sr);
            break;

        default: break;
    }
}


struct quip_in_t_
{
    quip_fmt_t fmt;
    union {
        quip_fastq_in_t* fastq;
        quip_sam_in_t*   sam;
        quip_quip_in_t*  quip;
    } x;
};


quip_in_t* quip_in_open(
              quip_reader_t reader,
              void*         reader_data,
              quip_fmt_t    fmt,
              quip_opt_t    opts)
{
    quip_in_t* in = malloc_or_die(sizeof(quip_in_t));
    in->fmt = fmt;

    switch (fmt) {
        case QUIP_FMT_FASTQ:
            in->x.fastq = quip_fastq_in_open(reader, reader_data, opts);
            break;

        case QUIP_FMT_BAM:
        case QUIP_FMT_SAM:
            in->x.sam = quip_sam_in_open(reader, reader_data, opts);
            break;

        case QUIP_FMT_QUIP:
            in->x.quip = quip_quip_in_open(reader, reader_data, opts);
            break;

        default: break;
    }

    return in;
}


void quip_in_close(quip_in_t* in)
{
    switch (in->fmt) {
        case QUIP_FMT_FASTQ:
            if (in->x.fastq) quip_fastq_in_close(in->x.fastq);
            break;

        case QUIP_FMT_BAM:
        case QUIP_FMT_SAM:
            if (in->x.sam) quip_sam_in_close(in->x.sam);
            break;

        case QUIP_FMT_QUIP:
            if (in->x.quip) quip_quip_in_close(in->x.quip);
            break;

        default: break;
    }

    free(in);
}


short_read_t* quip_read(quip_in_t* in)
{
    switch (in->fmt) {
        case QUIP_FMT_FASTQ:
            return quip_fastq_read(in->x.fastq);

        case QUIP_FMT_BAM:
        case QUIP_FMT_SAM:
            return quip_sam_read(in->x.sam);

        case QUIP_FMT_QUIP:
            return quip_quip_read(in->x.quip);

        default: break;
    }

    return NULL;
}


bool quip_pipe(quip_in_t* in, quip_out_t* out)
{
    short_read_t* sr = quip_read(in);
    if (sr == NULL) return false;
    else {
        quip_write(out, sr);
        return true;
    }
}
