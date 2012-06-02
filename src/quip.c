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
#include "seqmap.h"
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

static bool force_flag    = false;
static bool assembly_flag = false;
static bool stdout_flag   = false;

static enum {
    QUIP_CMD_CONVERT,
    QUIP_CMD_LIST,
} quip_cmd = QUIP_CMD_CONVERT;

static quip_fmt_t in_fmt  = QUIP_FMT_UNDEFINED;
static quip_fmt_t out_fmt = QUIP_FMT_UNDEFINED;

static bool force_in_fmt  = false;
static bool force_out_fmt = false;

static const char* ref_fn = NULL;

const char* fmt_suffix[] =
    {"", "", "fastq", "sam", "bam", "qp"};





static void print_help()
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
    quip_in_fname = fn;
    int fd = open(fn, O_RDONLY | O_NOCTTY | O_BINARY);

    if (fd == -1) {
        switch (errno) {
            case EACCES:
                quip_error("Permission denied.");
                break;

            default:
                quip_error("Error opening file.");
        }

        return NULL;
    }

    FILE* f = fdopen(fd, "rb");

    if (f == NULL) {
        quip_error("Error opening file.");
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
                    fprintf(stderr, "%s: %s: File already exists.\n", quip_prog_name, fn);
                    if (isatty(fileno(stdin))) {
                        fprintf(stderr, "Would you like to overwrite it (y or n)? ");
                        fflush(stderr);
                        overwrite = yesno();
                    }
                }

                if (overwrite) {
                    if (unlink(fn) == 0) return open_fout(fn);
                    fprintf(stderr, "%s: %s: Cannot overwrite file.\n", quip_prog_name, fn);
                }

                break;

            case EACCES:
                fprintf(stderr, "%s: %s: Permission denied.\n", quip_prog_name, fn);
                break;

            default:
                fprintf(stderr, "%s: %s: Error opening file.\n", quip_prog_name, fn);

        }

        return NULL;
    }

    FILE* f = fdopen(fd, "wb");

    if (f == NULL) {
        fprintf(stderr, "%s: %s: Error opening file.\n", quip_prog_name, fn);
        close(fd);
        return NULL;
    }

    return f;
}


/* Peek into a file to try to guess the type. */
static quip_fmt_t guess_file_format(const char* fn)
{
    char buf[1024];
    FILE* f = open_fin(fn);
    quip_fmt_t fmt;

    size_t n = fread(buf, sizeof(char), sizeof(buf) - 1, f);
    buf[n] = '\0';

    if (n == 0) fmt = QUIP_FMT_UNDEFINED;
    else {
        /* Check gzip header magic. */
        if (n >= 2 && memcmp(buf, "\037\213", 2) == 0) {
            fmt = QUIP_FMT_BAM;
        }
        /* Check quip header magic. */
        else if (n >= 6 && memcmp(buf, "\377QUIP\000", 6) == 0) {
            fmt = QUIP_FMT_QUIP;
        }
        else {
            /* Look at the leading character of consecutive lines
             * to choose between SAM and FASTQ. */
            char v, u = buf[0];
            char* line2 = strchr(buf, '\n');
            if (line2 == NULL) {
                fmt = QUIP_FMT_UNDEFINED;
            } 
            else {
                v = line2[1];

                if (u != '@' ||
                    (u == '@' && v == '@') ||
                    strchr(line2, '\t') != NULL) fmt = QUIP_FMT_SAM;
                else fmt = QUIP_FMT_FASTQ;
            }
       }
    }

    fclose(f);
    return fmt;
}



static int quip_cmd_convert(char** fns, size_t fn_count)
{
    if (stdout_flag) {
        SET_BINARY_MODE(stdout);
    }

    seqmap_t* ref = NULL;
    if (ref_fn != NULL) {
        ref = seqmap_alloc();
        seqmap_read_fasta(ref, ref_fn);
    }

    char* out_fn = NULL;

    FILE* fin;
    FILE* fout;

    quip_opt_t  opts;
    quip_in_t*  in;
    quip_out_t* out;
    quip_aux_t  aux;
    str_init(&aux.data);

    if (fn_count == 0) {
        quip_in_fname = "stdin";

        if (in_fmt == QUIP_FMT_UNDEFINED) {
            quip_warning("assuming input in FASTQ.");
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
            quip_error(
                "refusing to write compressed data to your terminal screen.\n\n"
                "Use -f is you really want to do this. (Hint: you don't.)");
        }

        in  = quip_in_open(block_reader,  (void*) stdin,  in_fmt, 0, ref);
        quip_get_aux(in, &aux);

        if (out_fmt == QUIP_FMT_QUIP && assembly_flag) opts = QUIP_OPT_QUIP_ASSEMBLY;
        else opts = 0;

        out = quip_out_open(block_writer, (void*) stdout, out_fmt, opts, &aux, ref);

        while (quip_pipe(in, out));

        quip_out_close(out);
        quip_in_close(in);

        fflush(stdout);
    }
    else {
        size_t i;
        for (i = 0; i < fn_count; ++i) {
            /* Determine the input format */
            if (in_fmt == QUIP_FMT_UNDEFINED) {
                in_fmt = guess_file_format(fns[i]);
                if (in_fmt == QUIP_FMT_UNDEFINED) {
                    quip_error("Unrecognized file format.");
                }
            }

            /* Determine the output format. */
            if (out_fmt == QUIP_FMT_UNDEFINED) {
                if (in_fmt == QUIP_FMT_FASTQ ||
                    in_fmt == QUIP_FMT_SAM ||
                    in_fmt == QUIP_FMT_BAM) {
                    out_fmt = QUIP_FMT_QUIP;
                }
                else {
                    /* Try guessing the output format from the file extension. */
                    char* b = strrchr(fns[i], '.');
                    char* a;
                    for (a = b - 1; a >= fns[i] && *a != '.'; --a);

                    if (a >= fns[i] && *a == '.') {
                        ++a;
                        if      (strcmp(a, "sam.qp") == 0) out_fmt = QUIP_FMT_SAM;
                        else if (strcmp(a, "bam.qp") == 0) out_fmt = QUIP_FMT_BAM;
                        else if (strcmp(a, "fastq.qp") == 0 ||
                                 strcmp(a, "fq.qp") == 0)  out_fmt = QUIP_FMT_FASTQ;
                    }

                    /* Our last resort guess of output format */
                    else if (ref == NULL) out_fmt = QUIP_FMT_FASTQ;
                    else             out_fmt = QUIP_FMT_SAM;
                }
            }

            fin = open_fin(fns[i]);
            in  = quip_in_open(block_reader, (void*) fin, in_fmt, 0, ref);

            quip_get_aux(in, &aux);

            /* Make a reasonable name for the output file. */
            if (stdout_flag) {
                fout = stdout;
            }
            else {
                if (out_fmt == QUIP_FMT_QUIP) {
                    asprintf(&out_fn, "%s.qp", fns[i]);
                }
                else if (in_fmt == QUIP_FMT_QUIP) {
                    size_t fn_len = strlen(fns[i]);

                    if (fn_len >= 3 && strcmp(fns[i] + fn_len - 3, ".qp") == 0) {
                        out_fn = malloc_or_die(fn_len - 3 + 1);
                        memcpy(out_fn, fns[i], fn_len - 3);
                        out_fn[fn_len] = '\0';
                    }
                    else {
                        asprintf(&out_fn, "%s.%s", fns[i], fmt_suffix[out_fmt]);
                    }
                }
                else {
                    size_t fn_len = strlen(fns[i]);
                    size_t in_suf_len = strlen(fmt_suffix[in_fmt]);

                    if (fn_len >= in_suf_len &&
                        strcmp(fns[i] + fn_len - in_suf_len, fmt_suffix[in_fmt]) == 0)
                    {
                        size_t out_suf_len = strlen(fmt_suffix[out_fmt]);
                        out_fn = malloc_or_die(fn_len - in_suf_len + 1 + out_suf_len);
                        memcpy(out_fn, fns[i], fn_len - in_suf_len);
                        memcpy(out_fn + fn_len - in_suf_len, fmt_suffix[out_fmt], out_suf_len);
                        out_fn[fn_len - in_suf_len + out_suf_len] = '\0';

                    }
                    else {
                        asprintf(&out_fn, "%s.%s", fns[i], fmt_suffix[out_fmt]);
                    }
                }

                fout = open_fout(out_fn);

                free(out_fn);
            }

            if (out_fmt == QUIP_FMT_QUIP && assembly_flag) opts = QUIP_OPT_QUIP_ASSEMBLY;
            else opts = 0;

            out = quip_out_open(block_writer, (void*) fout, out_fmt, opts, &aux, ref);

            while (quip_pipe(in, out));

            fflush(fout);

            quip_out_close(out);
            quip_in_close(in);

            if (fout != stdout) fclose(fout);
            fclose(fin);

            if (!force_in_fmt)  in_fmt  = QUIP_FMT_UNDEFINED;
            if (!force_out_fmt) out_fmt = QUIP_FMT_UNDEFINED;
        }
    }

    str_free(&aux.data);
    if (ref) seqmap_free(ref);

    return 0;
}


static void quip_print_list(const char* fn, quip_list_t* l)
{
    if (quip_verbose) {
        printf("%10"PRIu64"  "
               "%12"PRIu64"     "
               "%12"PRIu64"   "
               "%12"PRIu64"    "
               "%0.4f      "
               "%12"PRIu64"    "
               "%12"PRIu64"         "
               "%0.4f      "
               "%12"PRIu64"    "
               "%12"PRIu64"     "
               "%0.4f       "
               "%12"PRIu64"     "
               "%12"PRIu64"      "
               "%0.4f  "
               "%s"
               "\n",
               l->num_reads, l->num_bases,
               l->id_bytes[0],   l->id_bytes[1],
               (double) l->id_bytes[1] / (double) l->id_bytes[0],
               l->aux_bytes[0],  l->aux_bytes[1],
               (double) l->aux_bytes[1] / (double) l->aux_bytes[0],
               l->seq_bytes[0],  l->seq_bytes[1],
               (double) l->seq_bytes[1] / (double) l->seq_bytes[0],
               l->qual_bytes[0], l->qual_bytes[1],
               (double) l->qual_bytes[1] / (double) l->qual_bytes[0],
               fn);
    }
    else {
        uint64_t total_bytes[2];
        total_bytes[0] = l->id_bytes[0] + l->aux_bytes[0] + l->seq_bytes[0] + l->qual_bytes[0] + l->num_reads;
        total_bytes[1] = l->id_bytes[1] + l->aux_bytes[1] + l->seq_bytes[1] + l->qual_bytes[1] + l->header_bytes;

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
}


static int quip_cmd_list(char** fns, size_t fn_count)
{
    const char* fn;
    size_t fn_len;
    FILE* fin;
    size_t i;
    quip_list_t l;

    if (quip_verbose) {
        printf("     Reads         Bases  "
               "ID Uncompressed  ID Compressed  ID Ratio  "
               "Aux Uncompressed  Aux Compressed   Aux Ratio  "
               "Seq Uncompressed  Seq Compressed  Seq Ratio  "
               "Qual Uncompressed  Qual Compressed  Qual Ratio  "
               "Filename\n");
    }
    else {
        printf("     Reads         Bases  Uncompressed    Compressed   Ratio  Filename\n");
    }

    if (fn_count == 0) {
        quip_list(block_reader, stdin, &l);
        quip_print_list("stdin", &l);
    }
    else {
        for (i = 0; i < fn_count; ++i) {
            fn = fns[i];
            fn_len = strlen(fn);

            if (!force_flag && (fn_len < 3 || memcmp(fn + (fn_len - 3), ".qp", 3) != 0)) {
                quip_warning("unknown suffix -- ignored.");
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
        {"reference",  required_argument, NULL, 'r'},
        {"assembly-n", required_argument, NULL, 'n'},
        {"assembly",   no_argument      , NULL, 'a'},
        {"list",       no_argument, NULL, 'l'},
        {"test",       no_argument, NULL, 't'},
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
    quip_prog_name = argv[0];
    char* p;
    if ((p = strrchr(argv[0], '/')) != NULL) quip_prog_name = p + 1;
#if defined(WIN32) || defined(MSDOS)
    if ((p = strrchr(argv[0], '\\')) != NULL) quip_prog_name = p + 1;
#endif


    /* default to decompress when invoked under the name 'unquip' */
    if (strcmp(quip_prog_name, "unquip") == 0) {
        in_fmt = QUIP_FMT_QUIP;
    }
    else if (strcmp(quip_prog_name, "quipcat") == 0) {
        stdout_flag = true;
    }
    
    while (1) {
        opt = getopt_long(argc, argv, "i:o:r:n:ltacdfvhV", long_options, &opt_idx);

        if (opt == -1) break;

        switch (opt) {
            case 'i':
                in_fmt = parse_format(optarg);
                force_in_fmt = true;
                break;

            case 'o':
                out_fmt = parse_format(optarg);
                force_out_fmt = true;
                break;

            case 'r':
                ref_fn = optarg;
                break;

            case 'n':
                quip_assembly_n = strtoul(optarg, NULL, 10);
                assembly_flag = true;
                break;

            case 'l':
                quip_cmd = QUIP_CMD_LIST;
                break;

            case 't':
                in_fmt  = QUIP_FMT_QUIP;
                out_fmt = QUIP_FMT_NULL;
                break;

            case 'a':
                assembly_flag = true;
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

