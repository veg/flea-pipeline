#!/usr/bin/env python

"""Naive corrections to keep a sequence in-frame, relative to an
aligned in-frame reference.

The input file should contain pairs of aligned sequences and references.

Discards sequences that cannot be corrected, because they contain
insertions of length >3 that is not a multiple of three.

Usage:
  correct_shifts.py [options] <infile> <outfile>
  correct_shifts.py -h

Options:
  --keep            Do not discard sequences, even with bad inserts [default: False]
  --calns=<FILE>    File with run-length encoded alignment summaries.
  --discard=<FILE>  File to print discarded alignments.
  -v --verbose      Print summary [default: False]
  -h --help         Show this screen

"""

from itertools import groupby, repeat
from itertools import zip_longest
from functools import partial
import sys
import re

from docopt import docopt

from Bio import SeqIO
from Bio.SeqIO.FastaIO import FastaWriter

from util import grouper
from util import new_record_seq_str
from util import genlen
from util import write_to_handle
from util import nth


def first_index(target, it):
    """The first index in `it` where `target` appears; -1 if it is missing."""
    for idx, elt in enumerate(it):
        if elt == target:
            return idx
    return -1


def decode_caln(caln):
    """
    >>> decode_caln('3IDM2I2M')
    'IIIDMIIMM'

    """
    parts = re.findall(r'[0-9]*[IDM]', caln)
    if sum(map(len, parts)) != len(caln):
        raise Exception('badly formatted alignment encoding: {}'.format(caln))
    result = ''.join(part if len(part) == 1 else part[-1] * int(part[:-1])
                     for part in parts)
    return result


def alignment_slice(caln):
    """
    usearch trims leading and trailing insertions and deletions

    >>> alignment_slice("IMMMMMMMI")
    (2, 5)

    >>> alignment_slice("MMMDMMM")
    (0, 7)

    >>> alignment_slice("IMMMMM")
    (2, 5)

    >>> alignment_slice("MMMMMI")
    (0, 3)

    >>> alignment_slice("D9MD")
    (0, 9)

    >>> alignment_slice("D8MI")
    (0, 6)

    """
    aln = decode_caln(caln)
    if sum(1 for c in aln if c in 'MI') % 3 != 0:
        raise Exception('alignment has len(reference) not a multiple of 3')
    start_gaps = re.search('[M]', aln).start()
    end_gaps = re.search('[M]', aln[::-1]).start()

    start_insertions = sum(1 for c in aln[:start_gaps] if c == "I")
    end_insertions = sum(1 for c in aln[len(aln) - end_gaps:] if c == "I")

    # make aln match given alignment
    aln = aln[start_gaps:len(aln) - end_gaps]

    start_trim = (3 - start_insertions) % 3
    if start_trim > 0:
        # find start_trim'th match or insertion
        match = nth(re.finditer('[MI]', aln), start_trim)
        if match is None:
            raise Exception('deal with this')
        start = match.start()
    else:
        start = 0

    # find length of aligned target segment
    tlen = sum(1 for c in aln[start:] if c in 'MI')
    end_trim = (tlen - 3) % 3
    if end_trim > 0:
        # find end_trim'th match or insertion
        match = nth(re.finditer('[MI]', aln[::-1]), end_trim)
        if match is None:
            raise Exception('deal with this')
        stop = len(aln) - match.start()
    else:
        stop = len(aln)
    return start, stop


def correct_shifts(seq, ref, caln=None, gap_char=None, keep=False):
    """Correct frameshifts relative to a reference.

    keep: bool
        If True, keep sequences even if they could not be totally
        corrected. Otherwise, discard the whole sequence.

    """
    # FIXME: make this codon-aware
    if gap_char is None:
        gap_char = '-'
    start = 0
    if caln is not None:
        start, stop = alignment_slice(caln)
        seq = seq[start:stop]
        ref = ref[start:stop]
    if sum(1 for i in ref if i != gap_char) % 3 != 0:
        raise ValueError('len(reference) is not a multiple of 3')
    result = []
    index = start
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
                index += len(subseq)
                continue  # keep codon deletions
            if keep:
                result.append('N' * (len(subseq) % 3))
            else:
                return '', index
        else:  # insertion
            if len(subseq) % 3 == 0:
                result.append(subseq)  # keep codon insertions
            elif len(subseq) < 3:
                index += len(subseq)
                continue  # discard short insertions
            else:
                if keep:
                    # keep out-of-frame long insertions.
                    result.append(subseq)
                else:  # give up
                    return '', index
        index += len(subseq)
    if gap_char in result:
        raise Exception('gap character appeared in corrected result')
    result = ''.join(result)
    if not keep and len(result) % 3 != 0:
        raise Exception('shift correction failed')
    return result, -1


def correct_shifts_fasta(infile, outfile, calnfile=None,
                         discardfile=None,
                         alphabet=None, keep=True):
    """Correct all the pairs in a fasta file.

    Returns (n_seqs, n_fixed)

    """
    pairs = list(grouper(SeqIO.parse(infile, 'fasta', alphabet), 2))
    if calnfile is None:
        calns = repeat(None)
    else:
        with open(calnfile) as handle:
            calns = handle.read().strip().split()
    results = list(correct_shifts(seq.seq, ref.seq, caln, keep=keep)
                   for (seq, ref), caln in zip(pairs, calns))

    keep = list(new_record_seq_str(seq, result)
                for (result, _), (seq, ref) in zip(results, pairs) if result)
    SeqIO.write(keep, outfile, 'fasta')

    if discardfile is not None:
        discard = list(pair for (result, _), pair in zip(results, pairs)
                       if not result)
        for (result, index), (seq, ref) in zip(results, pairs):
            seq.description = "bad index: {}".format(index + 1)
            ref.description = "bad index: {}".format(index + 1)
        discard = list(s for pair in discard for s in pair)
        SeqIO.write(discard, discardfile, 'fasta')

    return len(results), len(keep)


def write_correction_result(n_seqs, n_fixed, handle):
    n_dropped = n_seqs - n_fixed
    percent = 100 * n_dropped / n_seqs
    do_close = False
    write_to_handle(handle, 'discarded {}/{} ({:.2f}%)'
                    ' sequences\n'.format(n_dropped, n_seqs, percent))


if __name__ == "__main__":
    args = docopt(__doc__)
    infile = args["<infile>"]
    outfile = args["<outfile>"]
    keep = args['--keep']
    calnfile = args['--calns']
    discardfile = args['--discard']
    n_seqs, n_fixed = correct_shifts_fasta(infile, outfile, discardfile=discardfile,
                                           calnfile=calnfile, keep=keep)
    if args['--verbose']:
        write_correction_result(n_seqs, n_fixed, sys.stdout)
