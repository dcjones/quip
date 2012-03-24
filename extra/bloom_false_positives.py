#!/usr/bin/env pypy

'''
Generate false-positive statistics for the d-left counting
bloom filter.
'''

import re
import subprocess
from random import randint
from sys    import stdout
# import argparse
# import numpy.random as npr



def fpr(xs, n, m):
    '''
    Given a collection of integers, measure the
    false positive rate of a n-by-m bloom filter.
    '''

    proc = subprocess.Popen(
        ['./bloom_stats', '-n', str(n), '-m', str(m)],
        stdout = subprocess.PIPE,
        stdin  = subprocess.PIPE)

    xs_data = '\n'.join(map(str, xs))
    (out_data, err_data) = proc.communicate(xs_data)

    mat = re.match(r'(\d+)\t(\d+)', out_data)
    assert mat

    u = int(mat.group(1))
    v = int(mat.group(2))
    return float(v - u) / v


def rand_int_set(n, low, high):
    '''
    Generate n random in [low, high]
    without replacement.
    '''

    xs = set()
    while len(xs) < n:
        xs.add(randint(low, high))

    return xs


def avg_fpr(c, n, m, low, high, rep = 10):
    '''
    Compute the false positive rate averaged
    over `rep` replicates.
    '''

    f = 0.0
    for _ in xrange(rep):
        xs = rand_int_set(c, low, high)
        f += fpr(xs, n, m)
    return f / rep


# k-mer size
k = 25

m = 8

# insertion counts
cs = [
        100,
       1000,
      10000,
     100000,
     200000,
     300000,
     400000,
     500000,
     600000,
     700000,
     800000,
     900000,
     910000,
     920000,
     930000,
     940000,
     950000,
     960000,
     970000,
     980000,
     990000,
    1000000]

# table sizes
ns = [1000000]


for c in cs:
    for n in ns:
        if c > n: continue
        print('{c}\t{n}\t{fpr}'.format(
            c = c, n = n, fpr = avg_fpr(c, n / (4 * m), m, 0, 4 ** 25 - 1)))
        stdout.flush()
