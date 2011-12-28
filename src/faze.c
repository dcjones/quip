
#include "assembler.h"
#include "kmer.h"
#include "misc.h"
#include "parse.h"
#include "qual.h"
#include <stdlib.h>

void qual_writer(void* param, uint8_t* data, size_t datalen)
{
    size_t* cnt = (size_t*) param;
    *cnt += datalen;
}


int main(int argc, char* argv[])
{
    kmer_init();

    FILE* fin;

    if (argc > 1) fin = fopen_or_die(argv[1], "r");
    else          fin = stdin;


    size_t k = 25;
    assembler_t* A = assembler_alloc(k);


    fastq_t* fqf = fastq_open(fin);
    seq_t* read = fastq_alloc_seq();

    size_t cnt = 0;

    size_t qual_bytes = 0;
    size_t qual_compressed_bytes = 0;
    qualenc_t* E = qualenc_alloc(qual_writer, &qual_compressed_bytes);


    while (fastq_next(fqf, read)) {
        if (++cnt % 100000 == 0) {
            fprintf(stdout, "%zu reads\n", cnt);
        }
        assembler_add_seq(A, read->seq.s, read->seq.n);
        qual_bytes += read->qual.n;
        /*qualenc_encode(E, read);*/
    }
    qualenc_finish(E);

    fastq_free_seq(read);
    fastq_close(fqf);

    qualenc_free(E);

    printf("%zu\t%zu\n", qual_bytes, qual_compressed_bytes);
    printf("%0.2f%% compression\n", 100.0 * (double) qual_compressed_bytes / (double) qual_bytes);


    assembler_write(A, stdout);
    assembler_free(A);

    kmer_free();

    return 0;
}

