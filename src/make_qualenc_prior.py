#!/usr/bin/env pypy

import argparse
import numpypy as np
from sys import stdout

max_cnt = 250

ap = argparse.ArgumentParser()
ap.add_argument('reads_fn', metavar = 'reads.fastq')
args = ap.parse_args()

xs = np.zeros(41 * 41, dtype = np.uint64)

f = open(args.reads_fn)
qq = 0
while True:
    id1  = f.readline()
    seq  = f.readline()
    id2  = f.readline()
    qual = f.readline()

    if id1 == '' or seq == '' or id2 == '' or qual == '':
        break

    for q in qual.strip():
        qq = (41 * qq + (ord(q) - 33)) % (41 * 41)
        xs[qq] += 1


for q in xrange(41):
    z = 0
    for u in xrange(41):
        z += xs[41 * q + u]

    d = 1 + z / max_cnt
    for u in xrange(41):
        xs[41 * q + u] = 1 + xs[41 * q + u] / d


stdout.write('#include "qualenc_prior.h"\n\n')
stdout.write('/* This file was generated automatically with "make_qualenc_prior.py".  */\n')
stdout.write('const uint16_t qualenc_prior[{0}] =\n'.format(41 * 41))
stdout.write('  { ')
for (i, x) in enumerate(xs[:-1]):
    if i > 0 and i % 11 == 0:
        stdout.write('\n    ')
    stdout.write(' {0:2d},'.format(int(x)))
stdout.write(' {0:2d} }};\n'.format(int(xs[-1])))


