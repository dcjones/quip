

#include <google/sparse_hash_map>
#include <cstdio>
#include <unistd.h>

typedef unsigned long kmer_t;

int main(int argc, char* argv[])
{
    size_t n = 1000000;
    int opt;

    while (true) {
        opt = getopt(argc, argv, "n:m:");

        if (opt == -1) break;

        switch (opt) {
            case 'n':
                n = (size_t) strtoul(optarg, NULL, 10);
                break;

            case '?':
                return EXIT_FAILURE;

            default:
                abort();
        }
    }

    google::sparse_hash_map<kmer_t, unsigned int> H(n);

    char line[512];

    kmer_t x;
    while (fgets(line, sizeof(line), stdin)) {
        x = (kmer_t) strtoul(line, NULL, 10);
        H[x] = 1;
    }

    return 0;
}
