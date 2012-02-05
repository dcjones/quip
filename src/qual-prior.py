#!/usr/bin/env pypy

'''
Generate a reasonable prior distribution over quality scores given some input
data.

Note that this closely mimics qualenc_encode/qualenc_decode. Don't change one of
these without changing the other.
'''

import argparse
import numpypy as np
from math import ceil
from sys import stdout, stderr

# number of pseudo-observations (i.e., the strength of the prior)
N = 500


pos_bins  = 8
qual_size = 41
qual_last = 40


ap = argparse.ArgumentParser()
ap.add_argument('reads_fn', metavar = 'reads.fastq')
args = ap.parse_args()

cnts = np.zeros(pos_bins * qual_size * qual_size * qual_size)

def tally(qs):
    n = len(qs)
    if n >= 1:
        cnts[qual_last * qual_size * qual_size +
             qual_last * qual_size +
             qs[0]] += 1

    if n >= 2:
        cnts[qual_last * qual_size * qual_size +
             qs[0] * qual_size +
             qs[1]] += 1

    if n >= 3:
        cnts[qs[0] * qual_size * qual_size +
             qs[1] * qual_size +
             qs[2]] += 1

    for i in xrange(3, n):
        cnts[(i * pos_bins) / (n + 1) * qual_size * qual_size * qual_size +
             max(qs[i - 3], qs[i - 2]) * qual_size * qual_size +
             qs[i - 1] * qual_size +
             qs[i]] += 1




f = open(args.reads_fn)
read_cnt = 0
while True:
    id1  = f.readline()
    seq  = f.readline()
    id2  = f.readline()
    qual = f.readline()

    if id1 == '' or seq == '' or id2 == '' or qual == '':
        break

    read_cnt += 1
    if read_cnt % 1000000 == 0:
        stderr.write('{0} reads ....\n'.format(read_cnt))

    qs = [ord(q) - 33 for q in qual.strip()]
    tally(qs)


# normalize each group
for i in xrange(pos_bins):
    for j in xrange(qual_size):
        for k in xrange(qual_size):
            idx0 = i * qual_size * qual_size * qual_size + \
                   j * qual_size * qual_size + \
                   k * qual_size

            z = 0
            for idx in xrange(idx0, idx0 + qual_size):
                z += cnts[idx]
            z = max(1, z)

            for idx in xrange(idx0, idx0 + qual_size):
                cnts[idx] = min(255, max(1, round(N * cnts[idx] / z)))


# serialize
stdout.write('#include "qual_prior.h"\n\n');
stdout.write('/* Automatically generated with qual-prior.py. Do not edit manually! */\n')
stdout.write('const uint8_t qual_prior[{0}] = \n'.format(len(cnts)))
stdout.write('    {')

for (i, c) in enumerate(cnts[:-1]):
    if i > 0 and i % 30 == 0:
        stdout.write('\n     ');
    stdout.write(' {0:3d},'.format(int(c)))

stdout.write(' {0:3d} }};\n'.format(int(cnts[-1])));


