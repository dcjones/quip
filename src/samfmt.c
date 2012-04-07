
#include "samfmt.h"
#include "misc.h"
#include "sam/sam.h"

struct quip_sam_out_t_
{
    samfile_t* f;
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

    return out;
}

void quip_sam_out_close(quip_sam_out_t* out)
{
    if (out) {
        samclose(out->f);
        free(out);
    }
}


struct quip_sam_in_t
{
    samfile_t* f;
};

