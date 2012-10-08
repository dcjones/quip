/*
 * This file is part of quip.
 *
 * Copyright (c) 2012 by Daniel C. Jones <dcjones@cs.washington.edu>
 *
 */

/* markov :
 * A sparse markov chain representation using the d-left counting bloom filter.
 */

#ifndef QUIP_MARKOV
#define QUIP_MARKOV

#include "kmer.h"
#include "ac.h"

typedef struct markov_t_ markov_t;


/* Initialize a new markov chain representaiton.
 *
 * Args:
 *   n: Allow up to this many k-mers in the sparse representation.
 *   k: Order of the sparse markov chain.
 *   k_catchall: Order of the fallback markov chain.
 *
 * Returns:
 *   An initialized markov_t pointer.
 */
markov_t* markov_create(size_t n, size_t k, size_t k_catchall);


/* Free a markov chain initialized with markov_create. */
void markov_free(markov_t* mc);


/* Encode a k-mer with the given markov chain, updating the model.
 *
 * Args:
 *   mc: A markov chain.
 *   ac: An arithmetic coder.
 *   ctx: The k-mer preceeding x, on which its probability is conditioned.
 *   x: A dinucleotide to encode.
 */
void markov_encode_and_update(markov_t* mc, ac_t* ac, kmer_t ctx, kmer_t x);


/* Decode a k-mer with a given markov chain, updating the model.
 *
 * Args:
 *   mc: A markov chain.
 *   ac: An arithmetic coder.
 *   ctx: The k-mer preceeding x, on which its probability is conditioned.
 *
 * Returns:
 *   A decoded dinucleotide.
 */
kmer_t markov_decode_and_update(markov_t* mc, ac_t* ac, kmer_t ctx);


#endif

