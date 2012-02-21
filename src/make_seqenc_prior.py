#!/usr/bin/env pypy
import numpypy as np

# Note: to run with normal python, replace the above two lines with:
# #!/usr/bin/env python
# import numpy as np

'''
Build a informed prior distribution over k-mers, given a genome sequenc.

This was used to generate ''seqenc_prior.c'. Regenerating it will create a
version of quip that is incompatible with other versions. So, it is not
recommended.
'''


import argparse
from sys  import stdout, stderr
from math import log

ap = argparse.ArgumentParser()
ap.add_argument('-k', default = 11, type = int)
ap.add_argument('genome_fn', metavar = 'genome.fa')
args = ap.parse_args()

max_cnt = 100

nucnum = {'A' : 0, 'C' : 1, 'G' : 2, 'T' : 3}


# Build complement tables
comp1 = [3, 2, 1, 0]
comp2 = [(comp1[x & 0x3] << 2)  | comp1[x >> 2] for x in xrange(0x10)]
comp4 = [(comp2[x & 0xf] << 4)  | comp2[x >> 4] for x in xrange(0x100)]
comp8 = [(comp4[x & 0xff] << 8) | comp4[x >> 8] for x in xrange(0x10000)]

def revcomp(x, k):
    y = (comp8[x & 0xffff] << 48) | \
        (comp8[(x >> 16) & 0xffff] << 32) | \
        (comp8[(x >> 32) & 0xffff] << 16) | \
        (comp8[(x >> 48) & 0xffff])

    return y >> (64 - (2 * k))



N = 4 ** (args.k + 2)
mask = N - 1
xs = np.zeros(N, dtype = np.uint32)

x = 0
for line in open(args.genome_fn):
    if line[0] == '>':
        stderr.write(line)
        continue
    vs = (nucnum[v] for v in line.strip().upper() if v != 'N')

    for v in vs:
        x = ((x << 2) | v) & mask
        xs[x] += 1
        xs[int(revcomp(x, args.k + 2))] += 1


for x in xrange(4 ** args.k):
    z = 0
    for u in xrange(16):
        z += xs[(x << 4) | u]

    # no count may exceed max_cnt
    q = 1 + z / max_cnt
    for u in xrange(16):
        xs[(x << 4) | u] = 1 + xs[(x << 4) | u] / q



# serialize
stdout.write('#include "seqenc_prior.h"\n\n')

stdout.write('/* This file was generated automatically with "make_seqenc_prior.py".  */\n')
stdout.write('const size_t seqenc_prior_k = {0};\n'.format(args.k));
stdout.write('const uint16_t seqenc_prior[{0}] =\n'.format(4 ** (args.k + 2)))
stdout.write('  { ')
for (i, f) in enumerate(xs[:-1]):
    if i > 0 and i % 16 == 0:
        stdout.write('\n    ')
    stdout.write(' {0:3d},'.format(int(f)))

stdout.write(' {0:3d} }};\n'.format(int(xs[-1])))

