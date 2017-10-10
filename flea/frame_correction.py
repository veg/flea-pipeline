#!/usr/bin/env python

"""Naive corrections to keep a sequence in-frame, relative to an
aligned in-frame reference.

The input file should contain pairs of aligned sequences and references.

Discards sequences that cannot be corrected, because they contain
insertions of length >3 that is not a multiple of three.

"""

from itertools import groupby, repeat
from itertools import zip_longest
from functools import partial
import re

import click

from Bio import SeqIO

from flea.util import grouper
from flea.util import new_record_seq_str
from flea.util import nth


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


def correct_frame(seq, ref, caln=None, gap_char=None, keep=True, deletion_strategy=None):
    """Correct frameshifts relative to a reference.

    keep: bool
        If True, keep sequences even if they could not be totally
        corrected. Otherwise, discard the whole sequence.

    """
    if deletion_strategy is None:
        deletion_strategy = 'reference'
    strategies = ['reference', 'n', 'discard']
    if deletion_strategy not in strategies:
        raise Exception('unknown strategy')

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
            elif len(subseq) == 1:
                if deletion_strategy == 'reference':
                    result.append(subref)
                elif deletion_strategy == 'n':
                    result.append('N')
                else:
                    # give up
                    return '', index
            elif keep:
                result.append('N' * (len(subseq) % 3))
            else:
                # give up
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


def correct_frames_fasta(infile, outfile, calnfile=None,
                         discardfile=None,
                         alphabet=None, keep=True, deletion_strategy=None):
    """Correct all the pairs in a fasta file.

    Returns (n_seqs, n_fixed)

    """
    pairs = list(grouper(SeqIO.parse(infile, 'fasta', alphabet), 2))
    if calnfile is None:
        calns = repeat(None)
    else:
        with open(calnfile) as handle:
            calns = handle.read().strip().split()
    results = list(correct_frame(seq.seq, ref.seq, caln,
                                 keep=keep, deletion_strategy=deletion_strategy)
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


def write_correction_result(n_seqs, n_fixed, outfile):
    n_dropped = n_seqs - n_fixed
    percent = 100 * n_dropped / n_seqs
    do_close = False
    with open(outfile, 'w') as handle:
        handle.write('discarded {}/{} ({:.2f}%)'
                     ' sequences\n'.format(n_dropped, n_seqs, percent))

@click.command()
@click.argument('infile')
@click.argument('outfile')
@click.option('--keep', help='eo not discard sequences, even with bad inserts')
@click.option('--deletion-strategy', help='correct single deletions')
@click.option('--calns', help='file with run-length encoded alignment summaries')
@click.option('--discard', help='file to print discarded alignments')
@click.option('--summary', help='file to print correction summaries')
def main(infile, outfile, keep, deletion_strategy, calns, discard, summary):
    n_seqs, n_fixed = correct_frames_fasta(infile, outfile, discardfile=discard,
                                           calnfile=calns, keep=keep,
                                           deletion_strategy=deletion_strategy)
    if summary:
        write_correction_result(n_seqs, n_fixed, summary)


if __name__ == "__main__":
    main()
