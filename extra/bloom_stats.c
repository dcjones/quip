
/*
  Empirical measurements of characteristics of
  the d-left counting bloom filter.

  The program will read integers from stdin and
  output the number of unique elements.
 */


#include "bloom.h"
#include <unistd.h>
#include <stdio.h>


int main(int argc, char* argv[])
{
    size_t n = 1000000;
    size_t m = 8;
    int opt;

    while (true) {
        opt = getopt(argc, argv, "n:m:");

        if (opt == -1) break;

        switch (opt) {
            case 'n':
                n = (size_t) strtoul(optarg, NULL, 10);
                break;

            case 'm':
                m = (size_t) strtoul(optarg, NULL, 10);
                break;

            case '?':
                return EXIT_FAILURE;

            default:
                abort();
        }
    }

    bloom_t* B = bloom_alloc(n, m);

    char line[512];

    kmer_t x;
    size_t count = 0, unique_count = 0;
    while (fgets(line, sizeof(line), stdin)) {
        x = (kmer_t) strtoul(line, NULL, 10);
        if (bloom_inc(B, x) == 1) ++unique_count;
        ++count;
    }

    bloom_free(B);

    printf("%zu\t%zu\n", unique_count, count);

    return 0;
}

