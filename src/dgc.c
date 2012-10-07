
/* Experiment with encoding against a de-bruijn graph. */

#include <stdio.h>

#include "quip.h"
#include "markov.h"


void byte_counter(void* param, const uint8_t* data, size_t datalen)
{
    size_t* count = (size_t*) param;
    *count += datalen;
}


int main(int argc, char* argv[])
{
    if (argc < 2) {
        fprintf(stderr, "Usage: dgc k < reads.fastq > data.out\n");
        return EXIT_FAILURE;
    }

    size_t k = atoi(argv[1]);

    size_t output_byte_count = 0;
    ac_t* ac = ac_alloc_encoder(byte_counter, &output_byte_count);

    markov_t* mc = markov_create(100000000, k);

    quip_in_t* in = quip_in_open_file(stdin, QUIP_FMT_FASTQ, 0, NULL);
    short_read_t* r;
    size_t input_byte_count = 0;
    size_t i;
    kmer_t x;
    kmer_t xmask = kmer_mask(k);
    while ((r = quip_read(in))) {
        input_byte_count += r->seq.n;
        x = 0;
        for (i = 0; i < r->seq.n; ++i) {
            x = ((x << 2) | chartokmer[r->seq.s[i]]) & xmask;
            markov_encode_and_update(mc, ac, x);
        }
    }
    quip_in_close(in);

    markov_free(mc);
    ac_flush_encoder(ac);
    ac_free(ac);

    printf("Input bytes: %zu\nOutput bytes: %zu\n",
           input_byte_count, output_byte_count);

    return EXIT_SUCCESS;
}


