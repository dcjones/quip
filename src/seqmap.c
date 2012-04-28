
#include "seqmap.h"
#include "quip.h"
#include "misc.h"
#include "crc64.h"
#include <stdlib.h>
#include <string.h>
#include <ctype.h>


typedef struct named_seq_t_
{
    char* seqname;
    twobit_t* seq;
    uint64_t  crc;
} named_seq_t;


static int named_seq_cmp(const void* a_, const void* b_)
{
    const named_seq_t* a = (named_seq_t*) a_;
    const named_seq_t* b = (named_seq_t*) b_;

    return strcmp(a->seqname, b->seqname);
}


struct seqmap_t_
{
    /* sequences */
    named_seq_t* seqs;

    /* number of elements allocated to seqs */
    size_t size;

    /* number of stored sequences */
    size_t n;
};


seqmap_t* seqmap_alloc()
{
    seqmap_t* M = malloc_or_die(sizeof(seqmap_t));
    M->seqs = NULL;
    M->size = 0;
    M->n    = 0;

    return M;
}


void seqmap_clear(seqmap_t* M)
{
    size_t i;
    for (i = 0; i < M->n; ++i) {
        twobit_free(M->seqs[i].seq);
        free(M->seqs[i].seqname);
    }

    M->n = 0;
}


void seqmap_free(seqmap_t* M)
{
    seqmap_clear(M);
    free(M);
}


static bool is_nt_char(char c)
{
    return c == 'a' || c == 'A' ||
           c == 'c' || c == 'C' ||
           c == 'g' || c == 'G' ||
           c == 't' || c == 'T' ||
           c == 'n' || c == 'N';
}


static void fasta_unexpected_char(char c)
{
    fprintf(stderr, "quip: error parsing fasta file: unexpected character '%c'\n", c);
    exit(EXIT_FAILURE);
}


static void check_unique(const seqmap_t* M)
{
    /* If there are more than one sequence going by the same name,
     * this will cause major problems. Check to make sure this is not
     * the case.
     */

    size_t i;
    for (i = 1; i < M->n; ++i) {
        if (strcmp(M->seqs[i].seqname, M->seqs[i - 1].seqname) == 0) {
            fprintf(stderr, "quip: reference contains multiple sequences of the same name: %s\n",
                            M->seqs[i].seqname);
            exit(EXIT_FAILURE);
        }
    }
}


void seqmap_read_fasta(seqmap_t* M, FILE* fasta_f)
{
    const size_t bufsize = 1024;
    char* buf = malloc_or_die(bufsize); buf[0] = '\0';
    char* next = buf;

    str_t seqname;
    str_init(&seqname);

    twobit_t* seq = NULL;

    /* The three parser states are:
         0 : reading seqname
         1 : reading sequence
         2 : reading sequence (line beginning)
    */
    int state = 2;

    while (true) {
        /* end of buffer */
        if (*next == '\0') {
            if (fgets(buf, bufsize, fasta_f) == NULL) {
                /* end of file */
                break;
            }
            else {
                next = buf;
                continue;
            }
        }

        else if (state == 0) {
            if (*next == '\n') {
                str_append_char(&seqname, '\0');
                if (quip_verbose) {
                    fprintf(stderr, "\treading %s...\n", seqname.s);
                }

                state = 2;

                if (seq != NULL) twobit_free_reserve(seq);
                seq = twobit_alloc();

                if (M->n >= M->size) {
                    M->size += 16;
                    M->seqs = realloc_or_die(M->seqs, M->size * sizeof(named_seq_t));
                }

                M->seqs[M->n].seq = seq;
                M->seqs[M->n].seqname = malloc_or_die(seqname.n);
                memcpy(M->seqs[M->n].seqname, seqname.s, seqname.n);
                M->n++;
            }
            else {
                str_append_char(&seqname, *next);
            }

            ++next;
        }
        else if (state == 1 || state == 2) {
            if (*next == '\n') {
                state = 2;
            }
            else if (is_nt_char(*next) && seq != NULL) {
                twobit_append_char(seq, *next);
                state = 1;
            }
            else if (state == 2 && *next == '>') {
                seqname.n = 0;
                seq = NULL;
                state = 0;
            }
            else {
                fasta_unexpected_char(*next);
            }

            ++next;
        }
    }

    if (seq != NULL) twobit_free_reserve(seq);

    qsort(M->seqs, M->n, sizeof(named_seq_t), named_seq_cmp);
    check_unique(M);

    str_free(&seqname);
    free(buf);
}


size_t seqmap_size(const seqmap_t* M)
{
    return M->n;
}


const twobit_t* seqmap_get(const seqmap_t* M, const char* seqname)
{
    if (M->n == 0) return NULL;

    /* sequences are stored sorted by name, so we lookup with binary search */
    size_t i = 0;
    size_t j = M->n - 1;
    size_t k;
    int c;

    while (i != j) {
        k = (i + j) / 2;
        c = strcmp(seqname, M->seqs[k].seqname);

        if (c < 0) i = k;
        else       j = k;
    }

    return strcmp(seqname, M->seqs[i].seqname) == 0 ? M->seqs[i].seq : NULL;
}


uint64_t seqmap_crc64(const seqmap_t* M)
{
    uint64_t crc = 0;
    size_t i;
    for (i = 0; i < M->n; ++i) {
        crc = crc64_update((uint8_t*) M->seqs[i].seqname, strlen(M->seqs[i].seqname), crc);
        crc = twobit_crc64_update(M->seqs[i].seq, crc);
    }

    return crc;
}

