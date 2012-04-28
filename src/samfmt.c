
#include "samfmt.h"
#include "misc.h"
#include "sam/sam.h"
#include "sam/kseq.h"
#include <inttypes.h>
#include <assert.h>

struct quip_sam_out_t_
{
    samfile_t* f;
    bam1_t* b;

    str_t sambuf;
    unsigned char* sambuf_next;
    tamFile sambuf_file;
};


size_t sambuf_reader(void* reader_data, uint8_t* data, size_t size)
{
    char** sambuf_next = (char**) reader_data;

    char* c = *sambuf_next;
    size_t i;
    for (i = 0; i < size && c[i]; ++i) {
        data[i] = (uint8_t) c[i];
    }

    *sambuf_next = c + i;

    return i;
}


quip_sam_out_t* quip_sam_out_open(
    quip_writer_t     writer,
    void*             writer_data,
    quip_opt_t        opts,
    const quip_aux_t* aux)
{
    quip_sam_out_t* out = malloc_or_die(sizeof(quip_sam_out_t));

    bool binary = (opts & QUIP_OPT_SAM_BAM) != 0;

    if (aux->fmt == QUIP_FMT_SAM || aux->fmt == QUIP_FMT_BAM) {
        out->f = samopen_out(writer, writer_data, binary, aux->aux);
    }
    else {
        bam_header_t* header = bam_header_init();
        const char* default_header_text = "@HD\tVN:1.0\tSO:unsorted\n";
        header->l_text = header->n_text = strlen(default_header_text);
        header->text = malloc_or_die(header->l_text + 1);
        memcpy(header->text, default_header_text, header->l_text + 1);

        out->f = samopen_out(writer, writer_data, binary, (void*) header);

        bam_header_destroy(header);
    }

    if (out->f == NULL) {
        fprintf(stderr, "Unable to open SAM/BAM output stream.\n");
        exit(EXIT_FAILURE);
    }

    out->b = bam_init1();
    str_init(&out->sambuf);
    out->sambuf_file = sam_open_in(sambuf_reader, (void*) &out->sambuf_next);

    return out;
}

void quip_sam_out_close(quip_sam_out_t* out)
{
    if (out) {
        samclose(out->f);
        bam_destroy1(out->b);
        sam_close(out->sambuf_file);
        str_free(&out->sambuf);
        free(out);
    }
}


void quip_sam_write(quip_sam_out_t* out, short_read_t* r)
{
    /* Packing SAM data into a bam1_t is complicated. Instead, I output a read
     * formatted in the SAM format to a temporary buffer and let samtools do the
     * conversion. */

    bool aligned = !(r->flags & BAM_FUNMAP);

    str_t* s = &out->sambuf; /* for convenience */
    s->n = 0;

    /* 1. qname */
    str_append(s, &r->id);

    /* 2. flag */
    str_reserve_extra(s, 12);
    s->n += snprintf((char*) s->s + s->n, 12, "\t%"PRIu32"\t", r->flags);

    /* 3. rname */
    if (aligned) {
        str_append(s, &r->seqname);
    }
    else {
        str_append_cstr(s, "*");
    }

    /* 4. pos */
    if (aligned) {
        str_reserve_extra(s, 12);
        s->n += snprintf((char*) s->s + s->n, 12, "\t%"PRIu32, r->pos + 1);
    }
    else {
        str_append_cstr(s, "\t0");
    }

    /* 5. mapq */
    str_reserve_extra(s, 6);
    s->n += snprintf((char*) s->s + s->n, 6, "\t%"PRIu8"\t",
            aligned ? r->map_qual : 255);

    /* 6. cigar */
    const char* cigar_chars = "MIDNSHP=X";
    size_t i;

    if (aligned && r->cigar.n > 0) {
        str_reserve_extra(s, r->cigar.n * 11 + 1);
        for (i = 0; i < r->cigar.n; ++i) {
            if (r->cigar.ops[i] > 8) continue;
            s->n += snprintf((char*) s->s + s->n, 12, "%"PRIu32"%c",
                        r->cigar.lens[i],
                        cigar_chars[r->cigar.ops[i]]);
        }
    }
    else {
        str_append_cstr(s, "*");
    }

    /* 7. rnext */
    if (aligned && r->mate_seqname.n > 0) {
        str_append_cstr(s, "\t");
        str_append(s, &r->mate_seqname);
    }
    else {
        str_append_cstr(s, "\t*");
    }

    /* 8. pnext */
    if (aligned) {
        str_reserve_extra(s, 12);
        s->n += snprintf((char*) s->s + s->n, 12, "\t%"PRIu32, r->mate_pos + 1);
    }
    else {
        str_append_cstr(s, "\t0");
    }

    /* 9. tlen */
    str_reserve_extra(s, 13);
    s->n += snprintf((char*) s->s + s->n, 13, "\t%"PRId32,
                aligned ? r->tlen : 0);

    /* 10. seq */
    if (r->seq.n > 0) {
        str_append_cstr(s, "\t");
        str_append(s, &r->seq);
    }
    else {
        str_append_cstr(s, "\t*");
    }

    /* 11. qual */
    if (r->qual.n > 0) {
        str_append_cstr(s, "\t");
        str_append(s, &r->qual);
    }
    else {
        str_append_cstr(s, "\t*");
    }

    /* 12. aux */
    uint8_t* a = (uint8_t*) r->aux.s;
    while (a < r->aux.s + r->aux.n) {
        uint8_t type, subtype, key[2];
        uint32_t array_len;
        key[0] = a[0]; key[1] = a[1];
        a += 2;
        type = *a;
        a += 1;

        str_reserve_extra(s, 5);
        s->n += snprintf((char*) s->s + s->n, 5, "\t%c%c:",
                (char) key[0], (char) key[1]);

        switch ((char) type) {
            case 'A':
                str_reserve_extra(s, 4);
                s->n += snprintf((char*) s->s + s->n, 4, "A:%c", (char) *a++);
                break;

            case 'C':
                str_reserve_extra(s, 6);
                s->n += snprintf((char*) s->s + s->n, 6, "i:%"PRIu8, *(uint8_t*) a);
                a++;
                break;

            case 'c':
                str_reserve_extra(s, 7);
                s->n += snprintf((char*) s->s + s->n, 7, "i:%"PRId8, *(int8_t*) a);
                a++;
                break;

            case 'S':
                str_reserve_extra(s, 8);
                s->n += snprintf((char*) s->s + s->n, 8, "i:%"PRIu16, *(uint16_t*) a);
                a += 2;
                break;

            case 's':
                str_reserve_extra(s, 9);
                s->n += snprintf((char*) s->s + s->n, 9, "i:%"PRId16, *(int16_t*) a);
                a += 2;
                break;

            case 'I':
                str_reserve_extra(s, 13);
                s->n += snprintf((char*) s->s + s->n, 13, "i:%"PRIu32, *(uint32_t*) a);
                a += 4;
                break;

            case 'i':
                str_reserve_extra(s, 14);
                s->n += snprintf((char*) s->s + s->n, 14, "i:%"PRId32, *(int32_t*) a);
                a += 4;
                break;

            case 'f':
                str_reserve_extra(s, 20);
                s->n += snprintf((char*) s->s + s->n, 20, "f:%g", *(float*) a);
                a += 4;
                break;

            case 'd':
                str_reserve_extra(s, 40);
                s->n += snprintf((char*) s->s + s->n, 40, "d:%lg", *(double*) a);
                a += 8;
                break;

            case 'Z':
            case 'H':
                str_reserve_extra(s, 3);
                s->n += snprintf((char*) s->s + s->n, 3, "%c:", type);
                while (*a) {
                    str_reserve_extra(s, 2);
                    s->s[s->n++] = *a++;
                }
                a++;
                break;

            case 'B':
                subtype = *a++;
                memcpy(&array_len, a, 4);
                a += 4;
                str_reserve_extra(s, 4);
                s->n += snprintf((char*) s->s + s->n, 4, "%c:%c", type, subtype);

                for (i = 0; i < array_len; ++i) {
                    str_append_cstr(s, ",");
                    switch (subtype) {
                        case 'c':
                            str_reserve_extra(s, 7);
                            s->n += snprintf((char*) s->s + s->n, 7, "i:%"PRId8, *(int8_t*) a);
                            a++;
                            break;

                        case 'C':
                            str_reserve_extra(s, 6);
                            s->n += snprintf((char*) s->s + s->n, 6, "i:%"PRIu8, *(uint8_t*) a);
                            a++;
                            break;

                        case 'S':
                            str_reserve_extra(s, 8);
                            s->n += snprintf((char*) s->s + s->n, 8, "i:%"PRIu16, *(uint16_t*) a);
                            a += 2;
                            break;

                        case 'i':
                            str_reserve_extra(s, 14);
                            s->n += snprintf((char*) s->s + s->n, 14, "i:%"PRId32, *(int32_t*) a);
                            a += 4;
                            break;

                        case 'I':
                            str_reserve_extra(s, 13);
                            s->n += snprintf((char*) s->s + s->n, 13, "i:%"PRIu32, *(uint32_t*) a);
                            a += 4;
                            break;

                        case 'f':
                            str_reserve_extra(s, 20);
                            s->n += snprintf((char*) s->s + s->n, 20, "f:%g", *(float*) a);
                            a += 4;
                            break;

                        default:
                            break;
                    }
                }

                break;

            default:
                break; // skip unknown aux type
        }
    }

    s->s[s->n] = '\0';

    out->sambuf_next = out->sambuf.s;
    sam_rewind(out->sambuf_file);
    sam_read1(out->sambuf_file, out->f->header, out->b);
    samwrite(out->f, out->b);
}


struct quip_sam_in_t_
{
    samfile_t* f;
    bam1_t* b;
    short_read_t r;
};


quip_sam_in_t* quip_sam_in_open(
                    quip_reader_t reader,
                    void*         reader_data,
                    ATTRIB_UNUSED quip_opt_t opts)
{
    quip_sam_in_t* in = malloc_or_die(sizeof(quip_sam_in_t));

    bool binary = (opts & QUIP_OPT_SAM_BAM) != 0;

    in->f = samopen_in(reader, reader_data, binary, NULL);

    if (in->f == NULL) {
        fprintf(stderr, "Unable to open SAM/BAM input stream.\n");
        exit(EXIT_FAILURE);
    }

    in->b = bam_init1();
    short_read_init(&in->r);

    return in;
}


void quip_sam_in_close(quip_sam_in_t* in)
{
    if (in) {
        samclose(in->f);
        bam_destroy1(in->b);
        short_read_free(&in->r);
        free(in);
    }
}


void quip_sam_get_aux(quip_sam_in_t* in, quip_aux_t* aux)
{
    aux->fmt = QUIP_FMT_SAM;
    aux->aux = (void*) in->f->header;
}


short_read_t* quip_sam_read(quip_sam_in_t* in)
{
    if (samread(in->f, in->b) <= 0) return NULL;

    str_copy_cstr(&in->r.id, (char*) bam1_qname(in->b), in->b->core.l_qname - 1);

    size_t readlen = (size_t) in->b->core.l_qseq;
    uint8_t* bamseq = bam1_seq(in->b);

    str_reserve(&in->r.seq, readlen + 1);
    size_t i;
    for (i = 0; i < readlen; ++i) {
        in->r.seq.s[i] = bam_nt16_rev_table[bam1_seqi(bamseq, i)];
    }
    in->r.seq.s[readlen] = '\0';
    in->r.seq.n = readlen;

    str_reserve(&in->r.qual, readlen + 1);
    uint8_t* bamqual = (uint8_t*) bam1_qual(in->b);
    for (i = 0; i < readlen; ++i) {
        in->r.qual.s[i] = bamqual[i] + 33;
    }
    in->r.qual.s[readlen] = '\0';
    in->r.qual.n = readlen;

    in->r.flags    = (uint32_t) in->b->core.flag;
    in->r.strand   = bam1_strand(in->b);
    if (in->r.strand) {
        /* BAM/SAM reorients reads to the the reference sequence */
        str_revcomp(in->r.seq.s, in->r.seq.n);
        str_rev(in->r.qual.s, in->r.qual.n);
    }

    in->r.pos      = in->b->core.pos;
    in->r.map_qual = in->b->core.qual;

    if (in->f->header && in->b->core.tid >= 0) {
        str_copy_cstr(
            &in->r.seqname,
            in->f->header->target_name[in->b->core.tid],
            strlen(in->f->header->target_name[in->b->core.tid]));
    }
    else in->r.seqname.n = 0;

    if (in->f->header && in->b->core.mtid >= 0) {
        if (in->b->core.mtid == in->b->core.tid) {
            in->r.mate_seqname.n = 0;
            str_append_cstr(&in->r.mate_seqname, "=");
        }
        else {
            str_copy_cstr(
                &in->r.mate_seqname,
                in->f->header->target_name[in->b->core.mtid],
                strlen(in->f->header->target_name[in->b->core.mtid]));
        }
    }
    else in->r.mate_seqname.n = 0;

    in->r.mate_pos = in->b->core.mpos;


    uint32_t* samcigar = bam1_cigar(in->b);
    size_t cigarlen = in->b->core.n_cigar;
    cigar_reserve(&in->r.cigar, cigarlen);

    for (i = 0; i < cigarlen; ++i) {
        in->r.cigar.ops[i]  = samcigar[i] & BAM_CIGAR_MASK;
        in->r.cigar.lens[i] = samcigar[i] >> BAM_CIGAR_SHIFT;
    }
    in->r.cigar.n = cigarlen;

    str_copy_cstr(&in->r.aux,
        (char*) bam1_aux(in->b),
        ((uint8_t*) in->b->data + in->b->data_len) - (uint8_t*) bam1_aux(in->b));


    return &in->r;
}


