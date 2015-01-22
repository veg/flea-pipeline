#!/usr/bin/env python

"""Naive corrections to keep a sequence in-frame, relative to an
aligned in-frame reference.

The input file should contain pairs of aligned sequences and references.

Discards sequences that cannot be corrected, because they contain
indels of length >3 that are not multiples of three.

Usage:
  correct_shifts.py [options] <infile> <outfile>
  correct_shifts.py -h

Options:
  -h --help              Show this screen

"""

import sys
from itertools import groupby
from itertools import zip_longest
from functools import partial

from docopt import docopt

from Bio import SeqIO
from Bio.SeqIO.FastaIO import FastaWriter

from util import grouper
from util import new_record_seq_str
from util import partition
from util import genlen


def first_index(target, it):
    """The first index in `it` where `target` appears; -1 if it is missing."""
    for idx, elt in enumerate(it):
        if elt == target:
            return idx
    return -1


def correct_shifts(seq, ref, gap_char=None):
    """Correct frameshifts relative to a reference."""
    if gap_char is None:
        gap_char = '-'
    result = []
    for k, g in groupby(zip_longest(seq, ref), key=partial(first_index, gap_char)):
        # k tells us where the gap is
        group = list(g)
        if any(c is None for subgroup in group for c in subgroup):
            raise ValueError('sequence and reference have different lengths')
        subseq, subref = tuple(''.join(chars) for chars in zip(*group))
        if k == -1:  # match
            result.extend(subseq)
        elif k == 0:  # deletion
            result.append('X' * (len(subseq) % 3))
        else:  # insertion
            if len(subseq) % 3 == 0:
                result.append(subseq)  # keep codon insertions
            elif len(subseq) < 3:
                continue  # discard other insertions
            else:
                return ''  # give up
    return ''.join(result)


def correct_shifts_fasta(infile, outfile):
    """Correct all the pairs in a fasta file.

    Returns (n_seqs, n_fixed)

    """
    ldict = {}
    pairs = grouper(SeqIO.parse(infile, 'fasta'), 2)
    results = (new_record_seq_str(seq, correct_shifts(seq.seq, ref.seq)) for seq, ref in pairs)

    results = genlen(results, ldict, 'n_seqs')
    fixed = genlen(filter(lambda r: r.seq, results), ldict, 'n_fixed')
    SeqIO.write(fixed, outfile, 'fasta')
    return ldict['n_seqs'], ldict['n_fixed']


if __name__ == "__main__":
    args = docopt(__doc__)
    infile = args["<infile>"]
    outfile = args["<outfile>"]
    n_seqs, n_fixed = correct_shifts_fasta(infile, outfile)
    n_dropped = n_seqs - n_fixed
    percent = 100 * n_dropped / n_seqs
    sys.stderr.write('correction done; discarded {}/{} ({:.2f}%)\n'.format(n_dropped, n_seqs, percent))
