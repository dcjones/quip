#!/usr/bin/env python

import argparse
import numpy as np
import Bio.SeqIO as SeqIO
from sys import stdout, stderr


ap = argparse.ArgumentParser()
ap.add_argument('-k', default = 4, type = int)
ap.add_argument('genome_fn', metavar = 'genome.fa')
args = ap.parse_args()


max_cnt = 16000

nucnum = {'A' : 0, 'C' : 1, 'G' : 2, 'T' : 3}


fs = np.zeros(4 ** (args.k + 2), dtype = int)
mask = (4 ** (args.k + 2)) - 1


def countkmers(xs):
    xs = (nucnum[u] for u in xs.upper() if u != 'N')
    x = 0
    for u in xs:
        x = ((x << 2) | u) & mask
        fs[x] += 1


for seqrec in SeqIO.parse(open(args.genome_fn), 'fasta'):
    print(seqrec.name)
    countkmers(seqrec.seq.reverse_complement().tostring())
    countkmers(seqrec.seq.tostring())



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


