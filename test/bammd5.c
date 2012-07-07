/*
 * This file is part of quip.
 *
 * Copyright (c) 2012 by Daniel C. Jones <dcjones@cs.washington.edu>
 *
 */

/*
  bammd5 : Compute MD5 hash on a SAM/BAM files.
  
  This is different than say `samtools view reads.bam | md5` because the SAM/BAM
  file is first parsed by quip. Additionally the order of optional data are
  ignored.
*/

#include "config.h"
#include "md5.h"
#include <quip.h>
#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>
#include <unistd.h>


void print_usage(FILE* out)
{
    fprintf(out,
            "bammd5 -- compute MD5 hashes on SAM and BAM files.\n\n"
            "Usage: cat reads.bam | bammd5 [-hsS]\n"
            "options:\n"
            "   -s, -S  input is SAM (by default BAM is assumed)\n" 
            "   -h      print this message\n");
}


int main(int argc, char* argv[])
{
    quip_fmt_t fmt = QUIP_FMT_BAM;

    int o;
    while ((o = getopt(argc, argv, "hsS")) != -1) {
        if (o == 's' || o == 'S') {
            fmt = QUIP_FMT_SAM;
        }
        else if (o == 'h') {
            print_usage(stdout);
            return 0;
        }
        else {
            print_usage(stderr);
            return 1;
        }
    }

    quip_in_t* in = quip_in_open_file(stdin, QUIP_FMT_BAM, 0, NULL);
    short_read_t* r;

    li_MD5_CTX md5ctx;
    li_MD5_Init(&md5ctx);
    uint32_t auxsize, i;

    samopt_t** opts = NULL;
    uint32_t opts_size = 0;

    while ((r = quip_read(in))) {
        li_MD5_Update(&md5ctx, r->id.s,              r->id.n);
        li_MD5_Update(&md5ctx, r->seq.s,             r->seq.n);
        li_MD5_Update(&md5ctx, r->qual.s,            r->qual.n);
        li_MD5_Update(&md5ctx, (void*) &r->flags,    sizeof(uint32_t));
        li_MD5_Update(&md5ctx, r->seqname.s,         r->seqname.n);
        li_MD5_Update(&md5ctx, (void*) &r->strand,   sizeof(uint8_t));
        li_MD5_Update(&md5ctx, (void*) &r->pos,      sizeof(uint32_t));
        li_MD5_Update(&md5ctx, (void*) &r->map_qual, sizeof(uint8_t));
        li_MD5_Update(&md5ctx, (void*) &r->cigar.n,  sizeof(uint32_t));
        li_MD5_Update(&md5ctx, r->cigar.ops,         r->cigar.n * sizeof(uint8_t));
        li_MD5_Update(&md5ctx, r->cigar.lens,        r->cigar.n * sizeof(uint32_t));
        li_MD5_Update(&md5ctx, r->mate_seqname.s,    r->mate_seqname.n);
        li_MD5_Update(&md5ctx, &r->mate_pos,         sizeof(uint32_t));
        li_MD5_Update(&md5ctx, &r->tlen,             sizeof(uint32_t));

        // hash optional fields without regard to order.
        auxsize = samopt_table_size(r->aux);
        li_MD5_Update(&md5ctx, &auxsize, sizeof(uint32_t));

        if (auxsize > opts_size) {
            opts_size = auxsize;
            opts = realloc(opts, opts_size * sizeof(samopt_t*));
        }

        samopt_table_dump_sorted(r->aux, opts);

        for (i = 0; i < auxsize; ++i) {
            li_MD5_Update(&md5ctx, &opts[i]->key,    2 * sizeof(uint8_t));
            li_MD5_Update(&md5ctx, &opts[i]->type,   sizeof(uint8_t));
            li_MD5_Update(&md5ctx, opts[i]->data->s, opts[i]->data->n);
        }
    }

    quip_in_close(in);
    free(opts);

    uint64_t md5[2];
    li_MD5_Final((unsigned char*) md5, &md5ctx);

    printf("%"PRIx64"%"PRIx64"\n", md5[0], md5[1]);

    return 0;
}