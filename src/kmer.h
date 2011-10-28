
#ifndef FAZER_KMER
#define FAZER_KMER

#include <stdint.h>

/* K-mers are encoded 2 bits per nucleotide in a 64 bit integer,
 * allowing up to 32.
 */
typedef uint64_t kmer_t;

/* nucleotide character to kmer */
kmer_t chartokmer(char);

/* nucleotide string to kmer */
kmer_t strtokmer(const char*);

/* complement */
kmer_t kmer_comp(kmer_t);

/* reverse complement */
kmer_t kmer_revcomp(kmer_t);

/* canonical */
kmer_t kmer_canonical(kmer_t);

#endif


