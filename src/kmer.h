
#ifndef FAZER_KMER
#define FAZER_KMER

#include <stdlib.h>
#include <stdint.h>

/* K-mers are encoded 2 bits per nucleotide in a 64 bit integer,
 * allowing up to k = 32.
 */
typedef uint64_t kmer_t;

/* nucleotide character to kmer */
kmer_t chartokmer(char);

/* kmer to nucleotide character */
char kmertochar(kmer_t);

/* nucleotide string to kmer */
kmer_t strtokmer(const char*);

/* kmer_t to character string */
void kmertostr(kmer_t, char*, size_t);

/* get a particular nucleotide from a sequence */
kmer_t kmer_get_nt(kmer_t* x, size_t i);

/* complement */
kmer_t kmer_comp(kmer_t, size_t k);

/* reverse complement */
kmer_t kmer_revcomp(kmer_t, size_t k);

/* canonical */
kmer_t kmer_canonical(kmer_t, size_t k);

/* hash functions for kmers */
uint64_t kmer_hash(kmer_t);
uint64_t kmer_hash_with_seed(kmer_t, uint64_t seed);


#endif


