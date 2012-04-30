#!/usr/bin/env pypy

"""
Generate run-time and memory usage statistics for the dlCBF and sparsehash.
"""

import numpypy as np
import re
import subprocess
import argparse
from random import randint


def timeit(cmd, input):
    proc = subprocess.Popen(
        ["gtime", "-f", "%e,%M"] + cmd,
        stdout = subprocess.PIPE,
        stderr = subprocess.PIPE,
        stdin  = subprocess.PIPE)

    (out, err) = proc.communicate(input)
    err = err.split("\n")
    print("\n".join(err[:-2]))
    print(out)

    mat = re.match(r"([\d\.]+),([\d\.]+)", err[-2])

    return (mat.group(1), mat.group(2))


def rand_ints(n, low, high):
    '''
    Generate n random in [low, high]
    with replacement.
    '''

    xs = np.empty(n, dtype=int)
    for i in xrange(n): xs[i] = randint(low, high)
    return xs


def rand_int_set(n, low, high):
    '''
    Generate n random in [low, high]
    without replacement.
    '''

    xs = set()
    while len(xs) < n:
        xs.add(randint(low, high))

    return xs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-n", default = 1000000, type = int)
    ap.add_argument("-k", default = 25, type = int)
    args = ap.parse_args()

    xs = rand_ints(args.n, 0, 4 ** args.k - 1)
    # xs = rand_int_set(args.n, 0, 4 ** args.k - 1)
    input = "\n".join(map(str, xs))

    N = int(0.7 * float(args.n / 8))
    (t,m) = timeit(["./bloom_stats", "-n", str(N)], input)
    print((t,m))

    (t,m) = timeit(["./sparsehash_stats", "-n", str(args.n)], input)
    print((t,m))


if __name__ == "__main__": main()
