
#ifndef FAZER_KMER
#define FAZER_KMER

#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>

/* K-mers are encoded 2 bits per nucleotide in a 64 bit integer,
 * allowing up to k = 32.
 */
typedef uint64_t kmer_t;


/* this function needs to be called when the program starts to build
 * reverse complement lookup tables, in particular. */
void kmer_init();
void kmer_free();


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

/* are fewer than all four nucleotides present in the kmer */
bool kmer_simple(kmer_t, size_t k);

/* hash functions for kmers */
uint64_t kmer_hash(kmer_t);
uint64_t kmer_hash_mix(uint64_t a, uint64_t b);

#endif


