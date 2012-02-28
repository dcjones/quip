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


static const char qual_first = 33;

const char* prog_name;

bool quick_flag      = false;
bool stdout_flag     = false;
bool decompress_flag = false;



void print_help()
{
    printf(
"Usage: quip [OPTION]... [FILE]...\n"
"Compress or decompress FASTQ sequence files with extreme prejudice.\n\n"
"  -d, --decompress   decompress\n"
"  -c, --stdout       write on standard output\n"
"  -q, --quick        compress much quicker, at a slightly lower compression ratio\n"
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

static size_t block_reader(void* param, uint8_t* data, size_t datalen)
{
    return fread(data, 1, datalen, (FILE*) param);
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

/* Note this function alters the input seq_t. Specifically, it
 * rescales the quality scores to printable characters. */
static void fastq_print(FILE* f, seq_t* seq)
{
    size_t i;
    for (i = 0; i < seq->qual.n; ++i) {
        seq->qual.s[i] += qual_first;
    }

    fprintf(f, "@%s\n%s\n+\n%s\n",
            seq->id1.s,
            seq->seq.s,
            seq->qual.s);
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
        C = quip_comp_alloc(block_writer, stdout, quick_flag);

        fq = fastq_open(stdin);

        while (fastq_next(fq, r)) {
            quip_comp_addseq(C, r);
        }

        fastq_close(fq);
        quip_comp_free(C);
    }
    else {
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
            }

            /* We do our own buffering, so we use unbuffered
             * streams to avoid unnecessary copying. */
            setvbuf(fin,  NULL, _IONBF, 0);
            setvbuf(fout, NULL, _IONBF, 0);

            C = quip_comp_alloc(block_writer, fout, quick_flag);

            fq = fastq_open(fin);

            while (fastq_next(fq, r)) {
                quip_comp_addseq(C, r);
            }

            quip_comp_flush(C);

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
    const char* fn;
    FILE* fin;
    FILE* fout;
    char* out_fn = NULL;
    size_t fn_len;
    size_t i;

    const size_t outbuf_size = 1048576;
    char* outbuf = malloc_or_die(outbuf_size * sizeof(char));

    seq_t* r = fastq_alloc_seq();
    quip_decompressor_t* D;

    if (fn_count == 0) {
        D = quip_decomp_alloc(block_reader, stdin);

        while (quip_decomp_read(D, r)) {
            fastq_print(stdout, r);
        }

        quip_decomp_free(D);

    }
    else {
        for (i = 0; i < fn_count; ++i) {
            fn = fns[i];
            fn_len = strlen(fn);

            if (fn_len < 3 || memcmp(fn + (fn_len - 3), ".qp", 3) != 0) {
                fprintf(stderr, "%s: %s: Unknown suffix -- ignored.\n",
                        prog_name, fn);
                continue;
            }

            fin = fopen_attempt(fn, "rb");
            if (!fin) continue;

            if (stdout_flag) fout = stdout;
            else {
                out_fn = malloc_or_die((fn_len - 2) * sizeof(char));
                memcpy(out_fn, fn, fn_len - 3);
                out_fn[fn_len - 2] = '\0';
                // TODO: test for clobber
                fout = fopen_attempt(out_fn, "wb");
                free(out_fn);
            }

            /* We do our own buffering for input, so we
             * use unbuffered streams to avoid unnecessary
             * copying. */
            setvbuf(fin,  NULL, _IONBF, 0);

            /* Output streams to benefit from buffering.
             * We use a very large to hopefully improve the
             * throughput of many small writes.
             */
            setvbuf(fin, outbuf, _IOFBF, outbuf_size);


            D = quip_decomp_alloc(block_reader, fin);

            while (quip_decomp_read(D, r)) {
                fastq_print(fout, r);
            }

            fclose(fin);

            quip_decomp_free(D);

            if (!stdout_flag) fclose(fout);

        }
    }


    fastq_free_seq(r);
    free(outbuf);

    return EXIT_SUCCESS;
}




int main(int argc, char* argv[])
{
    static struct option long_options[] = 
    {
        {"quick",      no_argument, NULL, 'q'},
        {"stdout",     no_argument, NULL, 'c'},
        {"decompress", no_argument, NULL, 'd'},
        {"verbose",    no_argument, NULL, 'v'},
        {"help",       no_argument, NULL, 'h'},
        {"version",    no_argument, NULL, 'V'},
        {NULL, 0, NULL, 0}
    };

    int opt, opt_idx;

    prog_name = argv[0];

    size_t prog_name_len = strlen(prog_name);

    /* default to decompress when invoked under the name 'unquip' */
    if (prog_name_len >= 6 && strcmp(argv[0] + (prog_name_len - 6), "unquip") == 0) {
        decompress_flag = true;
    }

    while (1) {
        opt = getopt_long(argc, argv, "qcdvhV", long_options, &opt_idx);

        if (opt == -1) break;

        switch (opt) {
            case 'q':
                quick_flag = true;
                break;

            case 'c':
                stdout_flag = true;
                break;

            case 'd':
                decompress_flag = true;
                break;

            case 'v':
                verbose = true;
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

