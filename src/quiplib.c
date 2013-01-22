
#include "config.h"
#include "quip.h"
#include "misc.h"
#include "seqmap.h"
#include "quipfmt.h"
#include "fastqfmt.h"
#include "samfmt.h"
#include "sam/bam.h"
#include <stdlib.h>
#include <stdarg.h>
#include <zlib.h>
#ifdef HAVE_LIBBZ2
#include <bzlib.h>
#endif

bool quip_verbose = false;
const char* quip_prog_name = "quip";
const char* quip_in_fname = "";


static void remove_output_file()
{
    /* TODO */
}


void quip_abort()
{
    remove_output_file();
    exit(EXIT_FAILURE);
}


void quip_error(const char* fmt, ...)
{
    va_list args;
    va_start(args, fmt);
    char* msg;
    vasprintf(&msg, fmt, args);
    va_end(args);

    fprintf(stderr, "\n%s: %s: %s\n",
            quip_prog_name, quip_in_fname, msg);

    free(msg);

    quip_abort();
}


void quip_warning(const char* fmt, ...)
{
    va_list args;
    va_start(args, fmt);
    char* msg;
    vasprintf(&msg, fmt, args);
    va_end(args);

    fprintf(stderr, "\n%s: %s: %s\n",
            quip_prog_name, quip_in_fname, msg);

    free(msg);
}


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


void str_reserve_extra(str_t* str, size_t size)
{
    if (str->n + size > str->size) {
        str->size = str->n + size;
        str->s = realloc_or_die(str->s, str->size * sizeof(char));
    }
}


void str_append(str_t* a, const str_t* b)
{
    str_reserve_extra(a, b->n);
    memcpy(a->s + a->n, b->s, b->n);
    a->n += b->n;
}


void str_append_cstr(str_t* a, const char* b)
{
    size_t len = strlen(b);
    str_reserve_extra(a, len);
    memcpy(a->s + a->n, b, len);
    a->n += len;
}


void str_append_char(str_t* a, char c)
{
    str_reserve_extra(a, 1);
    a->s[a->n++] = c;
}


void str_free(str_t* str)
{
    if (str) free(str->s);
}


void str_copy(str_t* dest, const str_t* src)
{
    str_reserve(dest, src->n + 1);
    if (src->s == NULL) dest->s[0] = '\0';
    else memcpy(dest->s, src->s, src->n);
    dest->n = src->n;
    dest->s[dest->n] = '\0';
}


void str_copy_cstr(str_t* dest, const char* src, size_t n)
{
    str_reserve(dest, n + 1);
    memcpy(dest->s, src, n);
    dest->s[n] = '\0';
    dest->n = n;
}


void str_memcpy(str_t* dest, const uint8_t* src, size_t n)
{
    str_reserve(dest, n);
    memcpy(dest->s, src, n);
    dest->n = n;
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


void cigar_copy(cigar_t* dest, const cigar_t* src)
{
    if (dest->size < src->size) {
        dest->size = src->size;
        dest->ops  = realloc_or_die(dest->ops,  dest->size * sizeof(uint8_t));
        dest->lens = realloc_or_die(dest->lens, dest->size * sizeof(uint32_t));
    }

    memcpy(dest->ops,  src->ops,  src->n * sizeof(uint8_t));
    memcpy(dest->lens, src->lens, src->n * sizeof(uint32_t));
    dest->n = src->n;
}


void short_read_init(short_read_t* sr)
{
    str_init(&sr->id);
    str_init(&sr->seq);
    str_init(&sr->qual);
    str_init(&sr->seqname);
    str_init(&sr->mate_seqname);
    sr->aux = samopt_table_alloc();
    cigar_init(&sr->cigar);
    sr->flags    = BAM_FUNMAP;
    sr->strand   = 0;
    sr->pos      = 0;
    sr->map_qual = 255;
    sr->mate_pos = 0;
    sr->tlen     = 0;
}


void short_read_free(short_read_t* sr)
{
    if (sr) {
        str_free(&sr->id);
        str_free(&sr->seq);
        str_free(&sr->qual);
        str_free(&sr->seqname);
        str_free(&sr->mate_seqname);
        samopt_table_free(sr->aux);
    }
}


void short_read_copy(short_read_t* dest, const short_read_t* src)
{
    str_copy(&dest->id,   &src->id);
    str_copy(&dest->seq,  &src->seq);
    str_copy(&dest->qual, &src->qual);
    str_copy(&dest->seqname, &src->seqname);
    str_copy(&dest->mate_seqname, &src->mate_seqname);
    samopt_table_copy(dest->aux, src->aux);
    cigar_copy(&dest->cigar, &src->cigar);
    dest->flags    = src->flags;
    dest->strand   = src->strand;
    dest->pos      = src->pos;
    dest->map_qual = src->map_qual;
    dest->mate_pos = src->mate_pos;
    dest->tlen     = src->tlen;
}



void quip_file_writer(void* param, const uint8_t* data, size_t datalen)
{
    fwrite(data, 1, datalen, (FILE*) param);
}


size_t quip_file_reader(void* param, uint8_t* data, size_t datalen)
{
    FILE* f = (FILE*) param;

    if (data == NULL) {
        if (fseek(f, datalen, SEEK_CUR) == 0) return datalen;
        else {
            /* The stream is not seekable, so we have
             * to read and discard datalen bytes. */
            const size_t bufsize = 4096;
            size_t readcnt = 0;
            uint8_t* buf = malloc_or_die(bufsize);
            size_t remaining = datalen;

            while (remaining > 0) {
                readcnt = fread(buf, 1, remaining < bufsize ? remaining : bufsize, f);
                if (readcnt == 0) break;
                else remaining -= readcnt;
            }

            free(buf);
            return datalen - remaining;
        }
    }
    else return fread(data, 1, datalen, (FILE*) param);
}


size_t quip_gzfile_reader(void* param, uint8_t* data, size_t datalen)
{
    gzFile gzf = (gzFile) param;

    if (data == NULL) {
        if (gzseek(gzf, datalen, SEEK_CUR) >= 0) {
            return datalen;
        }
        else {
            quip_error("Error reading gzipped file.");
            return 0; /* unreachable */
        }
    }
    else {
        int outlen = gzread(gzf, (void*) data, datalen);
        if (outlen < 0) {
            quip_error("Error reading gzipped file.");
        }
        return (size_t) outlen;
    }
}


#ifdef HAVE_LIBBZ2
size_t quip_bzfile_reader(void* param, uint8_t* data, size_t datalen)
{
    BZFILE* bzf = (BZFILE*) param;

    if (data == NULL) {
        quip_error("Seeking not supported by bzip2.");
        return 0; /* unreachable */
    }
    else {
        int outlen = BZ2_bzread(bzf, (void*) data, datalen);
        if (outlen < 0) {
            quip_error("Error reading bzip2 file.");
        }
        return (size_t) outlen;
    }
}
#endif


struct quip_out_t_
{
    quip_fmt_t fmt;
    union {
        quip_fastq_out_t* fastq;
        quip_sam_out_t*   sam;
        quip_quip_out_t*  quip;
        void*             null;
    } x;
};


quip_out_t* quip_out_open(
              quip_writer_t     writer,
              void*             writer_data,
              quip_fmt_t        fmt,
              quip_opt_t        opts,
              const quip_aux_t* aux,
              const seqmap_t*   ref)
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
            out->x.sam = quip_sam_out_open(writer, writer_data, opts, aux);
            break;

        case QUIP_FMT_QUIP:
            out->x.quip = quip_quip_out_open(writer, writer_data, opts, aux, ref);
            break;

        case QUIP_FMT_NULL:
            out->x.null = NULL;
            break;

        case QUIP_FMT_UNDEFINED:
            quip_error("Undefined format given.");

    }

    return out;
}


quip_out_t* quip_out_open_file(
              FILE*             file,
              quip_fmt_t        format,
              quip_opt_t        opts,
              const quip_aux_t* aux,
              const seqmap_t*   ref)
{
    return quip_out_open(quip_file_writer, (void*) file, format, opts, aux, ref);
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

        case QUIP_FMT_NULL:
            break;

        default:
            quip_error("Write called on an unsupported format.");
    }
}


struct quip_in_t_
{
    quip_fmt_t fmt;
    quip_filter_t filter;
    void* reader_data;

    union {
        quip_fastq_in_t* fastq;
        quip_sam_in_t*   sam;
        quip_quip_in_t*  quip;
    } x;
};


quip_in_t* quip_in_open(
              quip_reader_t   reader,
              void*           reader_data,
              quip_fmt_t      fmt,
              quip_filter_t   filter,
              quip_opt_t      opts,
              const seqmap_t* ref)
{
    quip_in_t* in = malloc_or_die(sizeof(quip_in_t));
    in->fmt = fmt;
    in->filter = filter;
    in->reader_data = reader_data;

    switch (fmt) {
        case QUIP_FMT_FASTQ:
            in->x.fastq = quip_fastq_in_open(reader, reader_data, opts);
            break;

        case QUIP_FMT_BAM:
            opts |= QUIP_OPT_SAM_BAM;

        case QUIP_FMT_SAM:
            in->x.sam = quip_sam_in_open(reader, reader_data, opts);
            break;

        case QUIP_FMT_QUIP:
            in->x.quip = quip_quip_in_open(reader, reader_data, opts, ref);
            break;

        default: break;
    }

    return in;
}


quip_in_t* quip_in_open_file(
              FILE*           file,
              quip_fmt_t      format,
              quip_filter_t   filter,
              quip_opt_t      opts,
              const seqmap_t* ref)
{
    if (filter == QUIP_FILTER_GZIP) {
        gzFile gzf = gzdopen(fileno(file), "rb");
        if (gzf == NULL) {
            quip_error("Cannot open gzipped file.");
        }
        return quip_in_open(quip_gzfile_reader, (void*) gzf, format,
                            filter, opts, ref);
    }
    else if (filter == QUIP_FILTER_BZIP2) {
#ifdef HAVE_LIBBZ2
        BZFILE* bzf = BZ2_bzdopen(fileno(file), "rb");
        if (bzf == NULL) {
            quip_error("Cannot open bzipped file.");
        }
        return quip_in_open(quip_bzfile_reader, (void*) bzf, format,
                            filter, opts, ref);
#else
        quip_error("Quip has not been compiled with bzip2 support.")
#endif
    }
    else if (filter == QUIP_FILTER_NONE) {
        return quip_in_open(quip_file_reader, (void*) file, format,
                            filter, opts, ref);
    }
    else {
        quip_error("Unsupported input filter.");
        return NULL; /* unreachable */
    }
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

    if (in->filter == QUIP_FILTER_GZIP) {
        gzclose((gzFile) in->reader_data);
    }
#ifdef HAVE_LIBBZ2
    else if (in->filter == QUIP_FILTER_BZIP2) {
        BZ2_bzclose((BZFILE*) in->reader_data);
    }
#endif

    free(in);
}


void quip_get_aux(quip_in_t* in, quip_aux_t* aux)
{
    aux->fmt = in->fmt;

    switch (in->fmt) {
        case QUIP_FMT_BAM:
        case QUIP_FMT_SAM:
            quip_sam_get_aux(in->x.sam, aux);
            break;

        case QUIP_FMT_QUIP:
            quip_quip_get_aux(in->x.quip, aux);
            break;

        default: break;
    }
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
