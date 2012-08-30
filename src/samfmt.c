
#include "samfmt.h"
#include "samopt.h"
#include "misc.h"
#include "sam/sam.h"
#include "sam/kseq.h"
#include <inttypes.h>
#include <assert.h>

struct quip_sam_out_t_
{
    samfile_t* f;
    bam1_t* b;
};

quip_sam_out_t* quip_sam_out_open(
    quip_writer_t     writer,
    void*             writer_data,
    quip_opt_t        opts,
    const quip_aux_t* aux)
{
    quip_sam_out_t* out = malloc_or_die(sizeof(quip_sam_out_t));

    bool binary = (opts & QUIP_OPT_SAM_BAM) != 0;

    bam_header_t* header = bam_header_init();

    const char* header_text;
    if (aux->fmt == QUIP_FMT_SAM || aux->fmt == QUIP_FMT_BAM) {
        header_text = (const char*) aux->data.s;
        header->n_text = header->l_text = aux->data.n;
    }
    else {
        header_text = "@HD\tVN:1.0\tSO:unsorted\n";
        header->n_text = header->l_text = header->n_text = strlen(header_text);
    }

    header->text = malloc_or_die(header->l_text + 1);
    memcpy(header->text, header_text, header->l_text);
    header->text[header->l_text] = '\0';
    sam_header_parse(header);

    out->f = samopen_out(writer, writer_data, binary, (void*) header);

    bam_header_destroy(header);

    if (out->f == NULL) {
        quip_error("Unable to open SAM/BAM output stream.");
    }

    out->b = bam_init1();

    return out;
}

void quip_sam_out_close(quip_sam_out_t* out)
{
    if (out) {
        samclose(out->f);
        bam_destroy1(out->b);
        free(out);
    }
}


static void copy_aux_field(str_t* s, const uint8_t** a_, uint8_t type)
{
    const uint8_t* a = *a_;
    uint8_t subtype;
    uint32_t array_len;

    switch ((char) type) {
        case 'A': case 'a': case 'c': case 'C':
            str_memcpy(s, a, 1);
            a += 1;
            break;

        case 'S': case 's':
            str_memcpy(s, a, 2);
            a += 2;
            break;

        case 'I': case 'i': case 'f':
            str_memcpy(s, a, 4);
            a += 4;
            break;

        case 'd':
            str_memcpy(s, a, 8);
            a += 8;
            break;

        case 'Z': case 'H':
            str_memcpy(s, a, strlen((char*) a) + 1);
            a += s->n;
            break;

        case 'B':
            subtype = *a;
            memcpy(&array_len, a + 1, 4);

            switch (subtype) {
                case 'A': case 'a': case 'c': case 'C':
                    str_memcpy(s, a, 5 + array_len);
                    break;

                case 'S': case 's':
                    str_memcpy(s, a, 5 + 2 * array_len);
                    break;

                case 'I': case 'i': case 'f':
                    str_memcpy(s, a, 5 + 4 * array_len);
                    break;

                case 'd':
                    str_memcpy(s, a, 5 + 8 * array_len);
                    break;

                default:
                    break;
            }
            a += s->n;
            break;

        default:
            break; // skip unknown aux type
    }

    *a_ = a;
}


static void read_aux(const uint8_t* src, size_t len, samopt_table_t* T)
{
    samopt_table_clear(T);

    samopt_t* opt;

    const uint8_t* a = src;
    while (a < src + len) {
        uint8_t type, key[2];
        key[0] = a[0]; key[1] = a[1];
        a += 2;
        type = *a;
        a += 1;

        opt = samopt_table_get(T, key);
        opt->type = type;

        copy_aux_field(opt->data, &a, type);
    }
}



static void bam_reserve_data(bam1_t* b, size_t size)
{
    if (b->m_data < (int) size) {
        b->m_data = (int) size;
        b->data = realloc_or_die(b->data, b->m_data);
    }
}


void quip_sam_write(quip_sam_out_t* out, short_read_t* r)
{
    /* for the sake of terseness */
    bam1_t*      b = out->b;
    bam1_core_t* c = &b->core;

    bool aligned      = !(r->flags & BAM_FUNMAP);
    bool mate_aligned = !(r->flags & BAM_FMUNMAP);

    size_t i, doff = 0;

    /* 1. qname */
    c->l_qname = r->id.n + 1;
    bam_reserve_data(b, c->l_qname);
    memcpy(b->data + doff, r->id.s, c->l_qname);
    doff += c->l_qname;

    /* 2. flags */
    c->flag = r->flags;

    /* 3. rname */
    c->tid = aligned ? bam_get_tid(out->f->header, (char*) r->seqname.s) : -1;

    /* 4. position */
    c->pos = aligned ? (int32_t) r->pos : -1;

    /* 5. mapq */
    c->qual = r->map_qual;

    /* 6. cigar (and bin) */
    bam_reserve_data(b, doff + 4 * r->cigar.n);
    c->n_cigar = r->cigar.n;
    for (i = 0; i < r->cigar.n; ++i) {
        bam1_cigar(b)[i] =
            (r->cigar.lens[i] << BAM_CIGAR_SHIFT) | r->cigar.ops[i];
    }
    doff += 4 * r->cigar.n;

    c->bin = bam_reg2bin(c->pos, bam_calend(c, bam1_cigar(b)));

    /* 7. rnext */
    c->mtid = (mate_aligned && r->mate_seqname.n > 0) ?
        bam_get_tid(out->f->header, (char*) r->mate_seqname.s) : -1;

    /* 8. pnext */
    c->mpos = mate_aligned ? (int32_t) r->mate_pos : -1;

    /* 9. tlen */
    c->isize = r->tlen;

    /* 10. seq and qual */
    c->l_qseq = r->qual.n;
    bam_reserve_data(b, doff + c->l_qseq + (c->l_qseq + 1) / 2);
    uint8_t* s = bam1_seq(b);
    memset(s, 0, (c->l_qseq + 1) / 2);

    /* Quip stores read in their original orientation. We need to reorient
     * for SAM/BAM. */
    if (r->strand) {
        uint8_t u;
        for (i = 0; i < r->seq.n; ++i) {
            u = complement(r->seq.s[r->seq.n - 1 - i]);
            s[i/2] |= bam_nt16_table[u] << (4 * (1 - i % 2));
        }
        doff += (c->l_qseq + 1) / 2;
        s = b->data + doff;

        if (r->qual.n == 0) {
            s[i] = 0xff;
            ++doff;
        }
        else {
            for (i = 0; i < r->seq.n; ++i) {
                u = r->qual.s[r->qual.n - 1 - i];
                s[i] = u - 33;
            }
            doff += c->l_qseq;
        }
    }
    else {
        for (i = 0; i < r->seq.n; ++i) {
            s[i/2] |= bam_nt16_table[(int) r->seq.s[i]] << (4 * (1 - i % 2));
        }
        doff += (c->l_qseq + 1) / 2;
        s = b->data + doff;

        if (r->qual.n == 0) {
            s[i] = 0xff;
            ++doff;
        }
        else {
            for (i = 0; i < r->seq.n; ++i) {
                s[i] = r->qual.s[i] - 33;
            }
            doff += c->l_qseq;
        }
    }

    /* 12. aux */
    samopt_table_bam_dump(r->aux, b);
    doff += b->l_aux;

    b->data_len = doff;

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
                    quip_opt_t    opts)
{
    UNUSED(opts);

    quip_sam_in_t* in = malloc_or_die(sizeof(quip_sam_in_t));

    bool binary = (opts & QUIP_OPT_SAM_BAM) != 0;

    in->f = samopen_in(reader, reader_data, binary, NULL);

    if (in->f == NULL) {
        quip_error("Unable to open SAM/BAM input stream.");
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
    /* SAM/BAM headers are stored in two parts. Sequences, with associated
     * lengths are stored in a hash table and remaining header fields are stored
     * in plaintext. But it's not quite that simple: the sequence ("@SQ") fields
     * may be repeated in the plaintext portion! Thus, samtools has graced us
     * with this ugly bit of code.
     */

    aux->fmt = QUIP_FMT_SAM;
    str_reserve(&aux->data, in->f->header->l_text);
    memcpy(aux->data.s, in->f->header->text, in->f->header->l_text);
    aux->data.n = in->f->header->l_text;

    bam_header_t* alt = bam_header_init();
    alt->l_text = in->f->header->l_text;
    alt->text   = in->f->header->text;
    sam_header_parse(alt);
    alt->l_text = 0;
    alt->text   = NULL;


    if (alt->n_targets) {
        if (alt->n_targets != in->f->header->n_targets) {
            quip_warning("Inconsistent number of target sequences in BAM header.");
        }
    }
    else {
        size_t max_bytes;
        int i;
        for (i = 0; i < in->f->header->n_targets; ++i) {
            max_bytes = 22 + strlen(in->f->header->target_name[i]);
            str_reserve_extra(&aux->data, max_bytes);
            aux->data.n += snprintf(
                    (char*) aux->data.s + aux->data.n,
                    max_bytes,
                    "@SQ\tSN:%s\tLN:%d\n",
                    in->f->header->target_name[i],
                    in->f->header->target_len[i]);
        }
    }

    bam_header_destroy(alt);
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

    read_aux(
        (const uint8_t*) bam1_aux(in->b),
        ((uint8_t*) in->b->data + in->b->data_len) - (uint8_t*) bam1_aux(in->b),
        in->r.aux);

    return &in->r;
}


