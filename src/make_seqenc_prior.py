#!/usr/bin/env pypy

import argparse
import numpypy as np
from sys import stdout, stderr


ap = argparse.ArgumentParser()
ap.add_argument('-k', default = 4, type = int)
ap.add_argument('genome_fn', metavar = 'genome.fa')
args = ap.parse_args()


max_cnt = 100

nucnum = {'A' : 0, 'C' : 1, 'G' : 2, 'T' : 3}


# a k-mer generator, given a fasta file
def kmers(k, f):
    mask = (4 ** k) - 1
    x = 0
    for line in f:
        line = line.strip()
        if len(line) == 0: continue
        if line[0] == '>':
            stderr.write('{0}\n'.format(line))
            continue

        for u in (nucnum[u] for u in line.upper() if u != 'N'):
            x = ((x << 2) | u) & mask
            yield x



fs = np.zeros(4 ** (args.k + 2), dtype = int)

for x in kmers(args.k + 2, open(args.genome_fn)):
    fs[x] += 1


for x in xrange(4 ** args.k):
    z = 0
    for u in xrange(16):
        z += fs[(x << 4) | u]

    # no count may exceed max_cnt
    q = 1 + z / max_cnt
    for u in xrange(16):
        fs[(x << 4) | u] = max(1, fs[(x << 4) | u] / q)


# serialize
stdout.write('#include "seqenc_prior.h"\n\n')

stdout.write('/* This file was generated automatically with "make_seqenc_prior.py".  */\n')
stdout.write('const size_t seqenc_prior_k = {0};\n'.format(args.k));
stdout.write('const uint16_t seqenc_prior[{0}] =\n'.format(4 ** (args.k + 2)))
stdout.write('  { ')
for (i, f) in enumerate(fs[:-1]):
    if i > 0 and i % 16 == 0:
        stdout.write('\n    ')
    stdout.write(' {0:3d},'.format(int(f)))

stdout.write(' {0:3d} }};\n'.format(int(fs[-1])))


