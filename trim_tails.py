#!/usr/bin/env python

"""
Trim poly-A tails from the front and end of sequences.

Sequences must be aligned, and all the same length'

Usage:
  trim_tails.py [options] <infile>
  trim_tails -h | --help

Options:
  -o --outfile=<STRING>  Name of output file
  -t --target=<CHAR>     Target letter [default: A]
  -n <INT>               Number of non-A's to terminate tail
  -h --help              Show this screen

"""

import re

import numpy as np
from Bio import SeqIO

from util import new_record_seq_str


def _rev(s):
    return "".join(reversed(s))


def _trim(seq, regex):
    match = regex.search(seq)
    return seq[match.start():]


def trim_tails(seq, target, n, regex=None):
    """

    >>> trim_tails("GAAAAAGGAA", "A", 2)
    'GG'

    >>> trim_tails("AAAGGAAA", "A", 2)
    'GG'

    """
    if regex is None:
        if len(target) != 1:
            raise Exception('"{}" is not a single character'.format(target))
        n = int(n)
        pattern = "((?!{}).){{{}}}".format(target, n)
        regex = re.compile(pattern)
    seq = seq.upper()
    result = _rev(_trim(_rev(seq), regex))
    return _trim(result, regex)


def trim_tails_file(infile, outfile, target, n):
    records = SeqIO.parse(infile, 'fasta')
    result = (new_record_seq_str(r, trim_tails(str(r.seq), target, n))
                                 for r in records)
    SeqIO.write(result, outfile, 'fasta')


if __name__ == "__main__":
    args = docopt(__doc__)
    filename = args["<infile>"]
    outfile = args["--outfile"]
    target = args["--target"]
    n = int(args["-n"])
    trim_tails_file(filename, outfile, target, n)

