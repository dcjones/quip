/*
 * This file is part of quip.
 *
 * Copyright (c) 2012 by Daniel C. Jones <dcjones@cs.washington.edu>
 *
 */


#include "config.h"
#include "quip.h"
#include "assembler.h"
#include "kmer.h"
#include "misc.h"
#include "parse.h"
#include <getopt.h>
#include <errno.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#if defined(MSDOS) || defined(OS2) || defined(WIN32) || defined(__CYGWIN__)
#  include <fcntl.h>
#  include <io.h>
#  define SET_BINARY_MODE(file) setmode(fileno(file), O_BINARY)
#else
#  define SET_BINARY_MODE(file)
#endif



const char* prog_name;

bool stdout_flag     = false;
bool decompress_flag = false;
bool keep_flag       = true;



void print_help()
{
    printf(
"Usage: quip [OPTION]... [FILE]...\n"
"Compress or decompress FASTQ sequence files with extreme prejudice.\n\n"
"  -c, --stdout       write on standard output, keep original files unchanged\n"
"  -d, --decompress   decompress\n"
"  -k, --keep         do not delete the input file(s)\n"
"  -v, --verbose      output lots of useless information\n"
"  -h, --help         print this message\n"
"  -V, --version      display program version\n\n"
"Report bugs to <dcjones@cs.washington.edu>.\n");
}


void print_version()
{
    printf("quip %s\n", VERSION);
}


static void block_writer(void* param, const uint8_t* data, size_t datalen)
{
    fwrite(data, 1, datalen, (FILE*) param);
}


/* Call fopen, and print an appropriate error message if need be. */
static FILE* fopen_attempt(const char* fn, const char* mode)
{
    FILE* f = fopen(fn, mode);
    if (f) return f;


    switch (errno) {
        case EACCES:
            fprintf(stderr, "%s: %s: Permission denied.\n", prog_name, fn);
            break;

        case ENOENT:
            fprintf(stderr, "%s: %s: No such file.\n", prog_name, fn);
            break;

        default:
            fprintf(stderr, "%s: %s: Error opening file.\n", prog_name, fn);
    }


    return NULL;
}


static int quip_compress(char** fns, size_t fn_count)
{
    if (stdout_flag) {
        SET_BINARY_MODE(stdout);
    }

    const char* fn;
    char* out_fn;
    FILE *fin, *fout;
    size_t i;

    seq_t* r = fastq_alloc_seq();
    fastq_t* fq;

    quip_compressor_t* C;

    if (fn_count == 0) {
        quip_write_header(stdout);
        C = quip_comp_alloc(block_writer, stdout);

        fq = fastq_open(stdin);

        while (fastq_next(fq, r)) {
            quip_comp_addseq(C, r);
        }

        fastq_close(fq);
        quip_comp_free(C);
    }
    else {
        if (stdout_flag) quip_write_header(stdout);

        for (i = 0; i < fn_count; ++i) {
            fn = fns[i];
            fin = fopen_attempt(fn, "rb");
            if (!fin) continue;

            if (stdout_flag) fout = stdout;
            else {
                out_fn = malloc_or_die((strlen(fn) + 4) * sizeof(char));
                sprintf(out_fn, "%s.qp", fn);
                fout = fopen_attempt(out_fn, "wb");
                free(out_fn);

                quip_write_header(fout);
            }

            C = quip_comp_alloc(block_writer, fout);

            fq = fastq_open(fin);

            while (fastq_next(fq, r)) {
                quip_comp_addseq(C, r);
            }


            fastq_close(fq);
            fclose(fin);

            quip_comp_free(C);
            if (!stdout_flag) fclose(fout);
        }
    }

    fastq_free_seq(r);

    return EXIT_SUCCESS;
}


static int quip_decompress(char** fns, size_t fn_count)
{
    // TODO
    return EXIT_SUCCESS;
}




int main(int argc, char* argv[])
{
    static struct option long_options[] = 
    {
        {"stdout",     no_argument, NULL, 'c'},
        {"decompress", no_argument, NULL, 'd'},
        {"keep",       no_argument, NULL, 'k'},
        {"help",       no_argument, NULL, 'h'},
        {"version",    no_argument, NULL, 'V'}
    };

    int opt, opt_idx;

    prog_name = argv[0];

    if (strcmp(argv[0], "unquip") == 0) decompress_flag = true;

    while (1) {
        opt = getopt_long(argc, argv, "cdkhV", long_options, &opt_idx);

        if (opt == -1) break;

        switch (opt) {
            case 'c':
                stdout_flag = true;
                break;

            case 'd':
                decompress_flag = true;
                break;

            case 'k':
                keep_flag = true;
                break;

            case 'h':
                print_help();
                return EXIT_SUCCESS;

            case 'V':
                print_version();
                return EXIT_SUCCESS;

            case '?':
                return EXIT_FAILURE;

            default:
                abort();
        }
    }

    /* initialize reverse complement lookup tables */
    kmer_init();

    int ret;
    if (decompress_flag) ret = quip_decompress(argv + optind, argc - optind);
    else                 ret = quip_compress(argv + optind, argc - optind);

    kmer_free();

    return ret;
}

