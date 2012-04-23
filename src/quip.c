/*
 * This file is part of quip.
 *
 * Copyright (c) 2012 by Daniel C. Jones <dcjones@cs.washington.edu>
 *
 */


#include "config.h"
#include "quip.h"
#include "quipfmt.h"
#include "fastqfmt.h"
#include "samfmt.h"
#include "assembler.h"
#include "kmer.h"
#include "misc.h"
#include <getopt.h>
#include <errno.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <inttypes.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/stat.h>
#ifndef O_BINARY
#define O_BINARY 0
#endif 

const char* prog_name;


bool force_flag  = false;
bool quick_flag  = false;
bool stdout_flag = false;

enum {
    QUIP_CMD_CONVERT,
    QUIP_CMD_LIST,
} quip_cmd = QUIP_CMD_CONVERT;

quip_fmt_t in_fmt  = QUIP_FMT_UNDEFINED;
quip_fmt_t out_fmt = QUIP_FMT_UNDEFINED;


void print_help()
{
    printf(
"Usage: quip [OPTION]... [FILE]...\n"
"Compress, decompress, or convert high-throughput\n"
"sequencing data with extreme prejudice.\n\n"
"Options:\n"
"  -i, --input=FORMAT, --from=FORMAT\n"
"                       input format (guessed by default)\n"
"  -o, --output=FORMAT, --to=FORMAT\n"
"                       output format (guessed by default)\n"
"  -d, --decompress     decompress (equivalent to '--input=quip')\n"
"  -t, --test           test compressed file integrity\n"
"  -l, --list           list total number of reads and bases\n"
"  -c, --stdout         write on standard output\n"
"  -f, --force          allow overwriting of output files, etc\n"
"  -q, --quick          compress quicker, at a lower compression ratio\n"
"  -v, --verbose        output a great deal of useless information\n"
"  -h, --help           print this message\n"
"  -V, --version        display program version\n\n"
"FORMAT is one of: quip, fastq, sam, bam\n"
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



static int quip_cmd_convert(char** fns, size_t fn_count)
{
    if (stdout_flag) {
        SET_BINARY_MODE(stdout);
    }

    quip_in_t*  in;
    quip_out_t* out;
    quip_aux_t  aux;


    if (fn_count == 0) {
        if (in_fmt == QUIP_FMT_UNDEFINED) {
            fprintf(stderr, "%s: stdin: Assuming input is FASTQ.\n", prog_name);
            in_fmt = QUIP_FMT_FASTQ;
        }

        if (out_fmt == QUIP_FMT_UNDEFINED) {
            if (in_fmt == QUIP_FMT_QUIP) {
                out_fmt = QUIP_FMT_SAM;
            }
            else out_fmt = QUIP_FMT_QUIP;
        }

        if (!force_flag && (out_fmt == QUIP_FMT_BAM || out_fmt == QUIP_FMT_QUIP) &&
            isatty(fileno(stdout)))
        {
            fprintf(stderr,
                "%s: refusing to write compressed data to your terminal screen.\n\n"
                "Use -f is you really want to do this. (Hint: you don't.)\n",
                prog_name);
            return EXIT_FAILURE;
        }

        in  = quip_in_open(block_reader,  (void*) stdin,  in_fmt, 0);
        quip_get_aux(in, &aux);
        out = quip_out_open(block_writer, (void*) stdout, out_fmt, 0, &aux);

        while (quip_pipe(in, out));

        quip_out_close(out);
        quip_in_close(in);

        fflush(stdout);
    }
    else {
        size_t i;
        for (i = 0; i < fn_count; ++i) {

        }
    }

    return 0;
}




#if 0
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

    quip_out_t* C;

    if (fn_count == 0) {
        C = quip_out_alloc(block_writer, stdout, quick_flag);

        fq = fastq_open(stdin);

        while (quip_out_readseq(C, fq));

        fastq_close(fq);
        quip_out_free(C);
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

            C = quip_out_alloc(block_writer, fout, quick_flag);

            fq = fastq_open(fin);

            while (quip_out_readseq(C, fq));

            quip_out_finish(C);

            fastq_close(fq);
            fclose(fin);

            quip_out_free(C);
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

    seq_t* r;
    quip_in_t* D;

    if (fn_count == 0) {
        D = quip_in_alloc(block_reader, stdin);

        while ((r = quip_in_read(D))) {
            if (!dryrun) fastq_print(stdout, r);
        }

        quip_in_free(D);
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


            D = quip_in_alloc(block_reader, fin);

            while ((r = quip_in_read(D))) {
                if (!dryrun) fastq_print(fout, r);
            }

            fclose(fin);

            quip_in_free(D);

            fflush(fout);
            if (!stdout_flag && !dryrun) fclose(fout);
        }
    }


    free(outbuf);

    return EXIT_SUCCESS;
}
#endif


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


static quip_fmt_t parse_format(const char* fmtstr)
{
    /* we might need something more sophisticated
     * when more formats are supported. */
    switch (tolower(fmtstr[0])) {
        case 'q': return QUIP_FMT_QUIP;
        case 'f': return QUIP_FMT_FASTQ;
        case 's': return QUIP_FMT_SAM;
        case 'b': return QUIP_FMT_BAM;
        default:  return QUIP_FMT_UNDEFINED;
    }
}


int main(int argc, char* argv[])
{
    static struct option long_options[] = 
    {
        {"input",      required_argument, NULL, 'i'},
        {"from",       required_argument, NULL, 'i'},
        {"output",     required_argument, NULL, 'o'},
        {"to",         required_argument, NULL, 'o'},
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
        in_fmt = QUIP_FMT_QUIP;
    }
    else if (strcmp(prog_name, "quipcat") == 0) {
        stdout_flag = true;
    }
    
    while (1) {
        opt = getopt_long(argc, argv, "ioltqcdfvhV", long_options, &opt_idx);

        if (opt == -1) break;

        switch (opt) {
            case 'i':
                in_fmt = parse_format(optarg);
                break;

            case 'o':
                out_fmt = parse_format(optarg);
                break;

            case 'l':
                quip_cmd = QUIP_CMD_LIST;
                break;

            case 't':
                in_fmt  = QUIP_FMT_QUIP;
                out_fmt = QUIP_FMT_NULL;
                break;

            case 'q':
                quick_flag = true;
                break;

            case 'c':
                stdout_flag = true;
                break;

            case 'd':
                in_fmt = QUIP_FMT_QUIP;
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
        case QUIP_CMD_CONVERT:
            ret = quip_cmd_convert(argv + optind, argc - optind);
            break;

       case QUIP_CMD_LIST:
            ret = quip_cmd_list(argv + optind, argc - optind);
            break;
    }

    kmer_free();

    return ret;
}

