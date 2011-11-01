
#include "assembler.h"
#include "misc.h"
#include "parse.h"
#include <stdlib.h>


int main(int argc, char* argv[])
{
    FILE* fin;

    if (argc > 1) fin = fopen_or_die(argv[1], "r");
    else          fin = stdin;


    size_t k = 25;
    assembler_t* A = assembler_alloc(k);


    fastq_t* fqf = fastq_open(fin);
    seq_t* read = fastq_alloc_seq();

    size_t cnt = 0;

    while (fastq_next(fqf, read)) {
        if (++cnt % 100000 == 0) {
            fprintf(stdout, "%zu reads\n", cnt);
        }
        assembler_add_seq(A, read->seq.s, read->seq.n);
    }

    fastq_free_seq(read);
    fastq_close(fqf);


    assembler_write(A, stdout);
    assembler_free(A);

    return 0;
}


