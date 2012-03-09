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
#include <stdio.h>
#include <string.h>
#include <inttypes.h>
#include <unistd.h>
#include <fcntl.h>
#ifndef O_BINARY
#define O_BINARY 0
#endif 

static const char qual_first = 33;

const char* prog_name;


bool force_flag      = false;
bool quick_flag      = false;
bool stdout_flag     = false;

enum {
    QUIP_CMD_COMPRESS,
    QUIP_CMD_DECOMPRESS,
    QUIP_CMD_LIST,
    QUIP_CMD_TEST
} quip_cmd = QUIP_CMD_COMPRESS;

void print_help()
{
    printf(
"Usage: quip [OPTION]... [FILE]...\n"
"Compress or decompress FASTQ sequence files with extreme prejudice.\n\n"
"  -d, --decompress   decompress\n"
"  -t, --test         test compressed file integrity\n"
"  -l, --list         list total number of reads and bases\n"
"  -c, --stdout       write on standard output\n"
"  -f, --force        allow overwriting of output files, etc\n"
"  -q, --quick        compress quicker, at a lower compression ratio\n"
"  -v, --verbose      output a great deal of useless information\n"
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

/* Prompt the user for a yes/no question. */
static bool yesno()
{
    int c = getchar();
    bool yes = c == 'y' || c == 'Y';
    while (c != '\n' && c != EOF) c = getchar();

    return yes;
}


/* Open an input file, or die trying */
static FILE* open_fin(const char* fn)
{
    int fd = open(fn, O_RDONLY | O_NOCTTY | O_BINARY);

    if (fd == -1) {
        switch (errno) {
            case EACCES:
                fprintf(stderr, "%s: %s: Permission denied.\n", prog_name, fn);
                break;

            default:
                fprintf(stderr, "%s: %s: Error opening file.\n", prog_name, fn);
        }

        return NULL;
    }

    FILE* f = fdopen(fd, "rb");

    if (f == NULL) {
        fprintf(stderr, "%s: %s: Error opening file.\n", prog_name, fn);
        close(fd);
        return NULL;
    }

    return f;
}


/* Open an output file, or die trying */
static FILE* open_fout(const char* fn)
{
    int fd = open(fn, O_WRONLY | O_CREAT | O_BINARY | O_EXCL,
                      S_IRUSR | S_IWUSR);
    bool overwrite = false;

    if (fd == -1) {
        switch (errno) {
            case EEXIST:
                if (force_flag) overwrite = true;
                else {
                    fprintf(stderr, "%s: %s: File already exists.\n", prog_name, fn);
                    if (isatty(fileno(stdin))) {
                        fprintf(stderr, "Would you like to overwrite it (y or n)? ");
                        fflush(stderr);
                        overwrite = yesno();
                    }
                }

                if (overwrite) {
                    if (unlink(fn) == 0) return open_fout(fn);
                    fprintf(stderr, "%s: %s: Cannot overwrite file.\n", prog_name, fn);
                }

                break;

            case EACCES:
                fprintf(stderr, "%s: %s: Permission denied.\n", prog_name, fn);
                break;

            default:
                fprintf(stderr, "%s: %s: Error opening file.\n", prog_name, fn);

        }

        return NULL;
    }

    FILE* f = fdopen(fd, "wb");

    if (f == NULL) {
        fprintf(stderr, "%s: %s: Error opening file.\n", prog_name, fn);
        close(fd);
        return NULL;
    }

    return f;
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


static int quip_cmd_compress(char** fns, size_t fn_count)
{
    if (stdout_flag) {
        SET_BINARY_MODE(stdout);
    }

    if (!force_flag && (stdout_flag || fn_count == 0) && isatty(fileno(stdout))) {
        fprintf(stderr,
            "%s: refusing to write compressed data to your terminal screen.\n\n"
            "Use -f is you really want to do this. (Hint: you don't.)\n",
            prog_name);
        return EXIT_FAILURE;
    }

    const char* fn;
    char* out_fn;
    FILE *fin, *fout;
    size_t i;

    fastq_t* fq;

    quip_compressor_t* C;

    if (fn_count == 0) {
        C = quip_comp_alloc(block_writer, stdout, quick_flag);

        fq = fastq_open(stdin);

        while (quip_comp_readseq(C, fq));

        fastq_close(fq);
        quip_comp_free(C);
    }
    else {
        for (i = 0; i < fn_count; ++i) {
            fn = fns[i];
            fin = open_fin(fn);
            if (!fin) continue;

            if (stdout_flag) fout = stdout;
            else {
                out_fn = malloc_or_die((strlen(fn) + 4) * sizeof(char));
                sprintf(out_fn, "%s.qp", fn);
                fout = open_fout(out_fn);
                free(out_fn);

                if (fout == NULL) {
                    fclose(fin);
                    continue;
                }
            }

            /* We do our own buffering, so we use unbuffered
             * streams to avoid unnecessary copying. */
            setvbuf(fin,  NULL, _IONBF, 0);
            setvbuf(fout, NULL, _IONBF, 0);

            C = quip_comp_alloc(block_writer, fout, quick_flag);

            fq = fastq_open(fin);

            while (quip_comp_readseq(C, fq));

            quip_comp_finish(C);

            fastq_close(fq);
            fclose(fin);

            quip_comp_free(C);
            if (!stdout_flag) fclose(fout);
        }
    }

    return EXIT_SUCCESS;
}


static int quip_cmd_decompress(char** fns, size_t fn_count, bool dryrun)
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
            if (!dryrun) fastq_print(stdout, r);
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

            fin = open_fin(fn);
            if (!fin) continue;

            if (dryrun) fout = NULL;
            else if (stdout_flag) fout = stdout;
            else {
                out_fn = malloc_or_die((fn_len - 2) * sizeof(char));
                memcpy(out_fn, fn, fn_len - 3);
                out_fn[fn_len - 2] = '\0';
                fout = open_fout(out_fn);
                free(out_fn);

                if (fout == NULL) {
                    fclose(fin);
                    continue;
                }
            }


            /* We do our own buffering for input, so we
             * use unbuffered streams to avoid unnecessary
             * copying. */
            setvbuf(fin,  NULL, _IONBF, 0);

            /* Output streams do benefit from buffering.
             * We use a very large to hopefully improve the
             * throughput of many small writes.
             */
            if (!dryrun) setvbuf(fout, outbuf, _IOFBF, outbuf_size);


            D = quip_decomp_alloc(block_reader, fin);

            while (quip_decomp_read(D, r)) {
                if(!dryrun) fastq_print(fout, r);
            }

            fclose(fin);

            quip_decomp_free(D);

            if (!stdout_flag && !dryrun) fclose(fout);

        }
    }


    fastq_free_seq(r);
    free(outbuf);

    return EXIT_SUCCESS;
}


static void quip_print_list(const char* fn, quip_list_t* l)
{
    const uint64_t read_overhead = 6;

    uint64_t total_bytes[2];
    total_bytes[0] = l->id_bytes[0] + l->seq_bytes[0] + l->qual_bytes[0] + read_overhead * l->num_reads;
    total_bytes[1] = l->id_bytes[1] + l->seq_bytes[1] + l->qual_bytes[1] + l->header_bytes;

    printf("%10"PRIu64"  "
           "%12"PRIu64"  "
           "%12"PRIu64"  "
           "%12"PRIu64"  "
           "%0.4f  "
           "%s"
           "\n",
           l->num_reads, l->num_bases,
           total_bytes[0], total_bytes[1],
           (double) total_bytes[1] / (double) total_bytes[0],
           fn);
}


static int quip_cmd_list(char** fns, size_t fn_count)
{
    const char* fn;
    size_t fn_len;
    FILE* fin;
    size_t i;
    quip_list_t l;

    printf("     Reads         Bases  Uncompressed    Compressed   Ratio  Filename\n");

    if (fn_count == 0) {
        quip_list(block_reader, stdin, &l);
        quip_print_list("stdin", &l);
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

            fin = open_fin(fn);
            if (!fin) continue;
            quip_list(block_reader, fin, &l);
            quip_print_list(fn, &l);
        }
    }

    return EXIT_SUCCESS;
}


int main(int argc, char* argv[])
{
    static struct option long_options[] = 
    {
        {"list",       no_argument, NULL, 'l'},
        {"test",       no_argument, NULL, 't'},
        {"quick",      no_argument, NULL, 'q'},
        {"stdout",     no_argument, NULL, 'c'},
        {"decompress", no_argument, NULL, 'd'},
        {"uncompress", no_argument, NULL, 'd'},
        {"force",      no_argument, NULL, 'f'},
        {"verbose",    no_argument, NULL, 'v'},
        {"help",       no_argument, NULL, 'h'},
        {"version",    no_argument, NULL, 'V'},
        {NULL, 0, NULL, 0}
    };

    int opt, opt_idx;

    /* determine the base program name */
    prog_name = argv[0];
    char* p;
    if ((p = strrchr(argv[0], '/')) != NULL) prog_name = p + 1;
#if defined(WIN32) || defined(MSDOS)
    if ((p = strrchr(argv[0], '\\')) != NULL) prog_name = p + 1;
#endif


    /* default to decompress when invoked under the name 'unquip' */
    if (strcmp(prog_name, "unquip") == 0) {
        quip_cmd = QUIP_CMD_DECOMPRESS;
    }
    else if (strcmp(prog_name, "quipcat") == 0) {
        quip_cmd = QUIP_CMD_DECOMPRESS;
        stdout_flag = true;
    }
    
    while (1) {
        opt = getopt_long(argc, argv, "ltqcdfvhV", long_options, &opt_idx);

        if (opt == -1) break;

        switch (opt) {
            case 'l':
                quip_cmd = QUIP_CMD_LIST;
                break;

            case 't':
                quip_cmd = QUIP_CMD_TEST;
                break;

            case 'q':
                quick_flag = true;
                break;

            case 'c':
                stdout_flag = true;
                break;

            case 'd':
                quip_cmd = QUIP_CMD_DECOMPRESS;
                break;

            case 'f':
                force_flag = true;
                break;

            case 'v':
                quip_verbose = true;
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
    switch (quip_cmd) {
        case QUIP_CMD_COMPRESS:
            ret = quip_cmd_compress(argv + optind, argc - optind);
            break;

        case QUIP_CMD_DECOMPRESS:
            ret = quip_cmd_decompress(argv + optind, argc - optind, false);
            break;

        case QUIP_CMD_LIST:
            ret = quip_cmd_list(argv + optind, argc - optind);
            break;

        case QUIP_CMD_TEST:
            ret = quip_cmd_decompress(argv + optind, argc - optind, true);
            break;
    }

    kmer_free();

    return ret;
}

