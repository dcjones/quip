
#include "samfmt.h"
#include "misc.h"
#include "sam/sam.h"

struct quip_sam_out_t_
{
    samfile_t* f;
    bam1_t* b;
};


quip_sam_out_t* quip_sam_out_open(
    quip_writer_t writer,
    void*         writer_data,
    quip_opt_t    opts)
{
    quip_sam_out_t* out = malloc_or_die(sizeof(quip_sam_out_t));

    bool binary = (opts & QUIP_OPT_SAM_BAM) != 0;

    /* TODO: we somehow need to retain the aux field. Otherwise
       we will lose the header after decompressing. */

    out->f = samopen_out(writer, writer_data, binary, NULL);

    if (out->f == NULL) {
        fprintf(stderr, "Unable to open SAM/BAM output stream.\n");
        exit(EXIT_FAILURE);
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


void quip_sam_write(quip_sam_out_t* out, short_read_t* r)
{
    // TODO
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
    }
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

    str_copy_cstr(&in->r.qual, (char*) bam1_qual(in->b),  in->b->core.l_qseq);

    in->r.flags  = (uint32_t) in->b->core.flag;
    in->r.strand = bam1_strand(in->b);
    in->r.pos    = in->b->core.pos;

    if (in->f->header && in->b->core.tid >= 0) {
        str_copy_cstr(
            &in->r.seqname,
            in->f->header->target_name[in->b->core.tid],
            strlen(in->f->header->target_name[in->b->core.tid]));
    }

    str_copy_cstr(&in->r.aux,
        (char*) bam1_aux(in->b),
        in->b->l_aux);


    /* TODO:
        cigar
        Other things? Look at the sam specification.
    */

    return &in->r;
}


