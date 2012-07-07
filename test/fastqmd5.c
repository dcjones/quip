/*
 * This file is part of quip.
 *
 * Copyright (c) 2012 by Daniel C. Jones <dcjones@cs.washington.edu>
 *
 */

/*
  fastqmd5 : Compute MD5 hash on a fastq file.
  
  This is different than say `cat reads.fastq | md5` because the fastq
  is first parsed by quip. This exclude the second id (following the '@'
  symbol) present in some fastq files.
*/


#include "config.h"
#include "md5.h"
#include <quip.h>
#include <stdio.h>
#include <inttypes.h>
#include <unistd.h>


void print_usage(FILE* out)
{
    fprintf(out,
            "fastqmd5 -- compute MD5 hashes on FASTQ files.\n\n"
            "Usage: cat reads.fastq | fastqmd5 [-h]\n"
            "options:\n"
            "   -h      print this message\n");
}


int main(int argc, char* argv[])
{
    int o;
    while ((o = getopt(argc, argv, "hsS")) != -1) {
        if (o == 'h') {
            print_usage(stdout);
            return 0;
        }
        else {
            print_usage(stderr);
            return 1;
        }
    }


    quip_in_t* in = quip_in_open_file(stdin, QUIP_FMT_FASTQ, 0, NULL);
    short_read_t* r;

    li_MD5_CTX md5ctx;
    li_MD5_Init(&md5ctx);

    while ((r = quip_read(in))) {
        li_MD5_Update(&md5ctx, r->id.s,   r->id.n);
        li_MD5_Update(&md5ctx, r->seq.s,  r->seq.n);
        li_MD5_Update(&md5ctx, r->qual.s, r->qual.n);
    }

    uint64_t md5[2];
    li_MD5_Final((unsigned char*) md5, &md5ctx);

    printf("%016"PRIx64"%016"PRIx64"\n", md5[0], md5[1]);

    quip_in_close(in);

    return 0;
}
