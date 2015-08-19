#!/usr/bin/env python

"""Naive corrections to keep a sequence in-frame, relative to an
aligned in-frame reference.

The input file should contain pairs of aligned sequences and references.

Discards sequences that cannot be corrected, because they contain
insertions of length >3 that is not a multiples of three.

Usage:
  correct_shifts.py [options] <infile> <outfile>
  correct_shifts.py -h

Options:
  --keep                 Do not discard sequences, even with bad inserts [default: False]
  -v --verbose           Print summary [default: False]
  -h --help              Show this screen

"""

from itertools import groupby, repeat
from itertools import zip_longest
from functools import partial
import sys

from docopt import docopt

from Bio import SeqIO
from Bio.SeqIO.FastaIO import FastaWriter

from util import grouper
from util import new_record_seq_str
from util import partition
from util import genlen
from util import write_to_handle


def first_index(target, it):
    """The first index in `it` where `target` appears; -1 if it is missing."""
    for idx, elt in enumerate(it):
        if elt == target:
            return idx
    return -1


def nth_nongap(seq, n, gap_char=None, reverse=False):
    """

    >>> nth_nongap('--a--bcd', 1)
    2

    >>> nth_nongap('--abc---', 1, reverse=True)
    4

    >>> nth_nongap('a--bcd', 1)
    0

    >>> nth_nongap('a--bcd', 2)
    3

    >>> nth_nongap('a--bc-d', 1, reverse=True)
    6

    >>> nth_nongap('a--bc-d', 2, reverse=True)
    4

    """
    if gap_char is None:
        gap_char = '-'
    if reverse:
        seq = list(reversed(seq))
    count = 0
    found = False
    result = -1
    for i, char in enumerate(seq):
        if char != gap_char:
            count += 1
        if count == n:
            found = True
            result = i
            break
    if not found:
        raise Exception('n is too large')
    if reverse:
        result = len(seq) - result - 1
    return result


def alignment_slice(ref, offsets=None, gap_char=None):
    """
    offsets are [start, stop]

    >>> alignment_slice("aaccc---ggg", [1, 8])
    (2, 11)

    >>> alignment_slice("aaa--bbb", [0, 5])
    (0, 8)

    >>> alignment_slice("aaa--bb", [0, 4])
    (0, 3)

    >>> alignment_slice("aabbbcccd", [1, 9])
    (2, 8)

    >>> alignment_slice("abbbcccdd", [2, 10])
    (1, 7)

    >>> alignment_slice("abbbcccddd", [2, 11])
    (1, 10)

    >>> alignment_slice("---aaa---", [0, 2])
    (0, 9)

    """
    # FIXME: leading and trailing insertions relative to reference
    if gap_char is None:
        gap_char = '-'
    start = 0
    stop = len(ref)
    if offsets is None:
        return start, stop
    ungapped_len = sum(1 for c in ref if c != gap_char)
    if offsets[1] - offsets[0] + 1 != ungapped_len:
        raise Exception('offsets are wrong')
    trim_start = (3 - (offsets[0] % 3)) % 3
    start = nth_nongap(ref, trim_start + 1, gap_char)

    trim_stop = (offsets[1] + 1) % 3
    stop = nth_nongap(ref, trim_stop + 1, gap_char, reverse=True) + 1

    new_len = len(list(c for c in ref[start:stop] if c != gap_char))
    if (new_len % 3) != 0:
        print(len(ref), ref)
        print(offsets)
        print(start, stop)
        print(sum(1 for c in ref[start:stop] if c != gap_char))
        raise Exception('alignment_slice failed to make an in-frame reference')
    return start, stop


def correct_shifts(seq, ref, offsets=None, gap_char=None, keep=False):
    """Correct frameshifts relative to a reference.

    keep: bool
        If True, keep sequences even if they could not be totally
        corrected. Otherwise, discard the whole sequence.

    """
    # FIXME: make this codon-aware
    if gap_char is None:
        gap_char = '-'
    start, stop = alignment_slice(ref, offsets, gap_char)
    seq = seq[start:stop]
    ref = ref[start:stop]
    if sum(1 for c in ref if c != gap_char) % 3 != 0:
        raise Exception('len(reference) is not a multiple of 3')
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
            if len(subseq) % 3 == 0:
                continue  # keep codon deletions
            if keep:
                result.append('X' * (len(subseq) % 3))
            else:
                return ''
        else:  # insertion
            if len(subseq) % 3 == 0:
                result.append(subseq)  # keep codon insertions
            elif len(subseq) < 3:
                continue  # discard short insertions
            else:
                if keep:
                    # keep out-of-frame long insertions.
                    result.append(subseq)
                else:  # give up
                    return ''
    if gap_char in result:
        raise Exception('gap character appeared in corrected result')
    if len(result) % 3 != 0:
        raise Exception('shift correction failed')
    return ''.join(result)


def correct_shifts_fasta(infile, outfile, offsetfile=None, alphabet=None, keep=True):
    """Correct all the pairs in a fasta file.

    Returns (n_seqs, n_fixed)

    """
    ldict = {}
    pairs = grouper(SeqIO.parse(infile, 'fasta', alphabet), 2)
    if offsetfile is None:
        offsets = repeat(None)
    else:
        with open(offsetfile) as handle:
            offsets = grouper(map(int, handle.read().split()), 2)
    results = (new_record_seq_str(seq, correct_shifts(seq.seq, ref.seq, off, keep=keep))
               for (seq, ref), off in zip(pairs, offsets))

    results = genlen(results, ldict, 'n_seqs')
    fixed = genlen(filter(lambda r: r.seq, results), ldict, 'n_fixed')
    SeqIO.write(fixed, outfile, 'fasta')
    return ldict['n_seqs'], ldict['n_fixed']


def write_correction_result(n_seqs, n_fixed, handle):
    n_dropped = n_seqs - n_fixed
    percent = 100 * n_dropped / n_seqs
    do_close = False
    write_to_handle(handle, 'discarded {}/{} ({:.2f}%)'
                    ' sequences'.format(n_dropped, n_seqs, percent))


if __name__ == "__main__":
    args = docopt(__doc__)
    infile = args["<infile>"]
    outfile = args["<outfile>"]
    keep = args['<--keep-inserts>']
    n_seqs, n_fixed = correct_shifts_fasta(infile, outfile, keep)
    if args['--verbose']:
        write_correction_result(n_seqs, n_fixed, sys.stdout)
