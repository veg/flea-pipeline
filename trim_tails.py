#!/usr/bin/env python

"""
Trim poly-A tails from the start or end of sequences.

Uses a simple HMM with unweighted transitions.

Usage:
  trim_tails.py [options] <head_p> <tail_p> <infile>
  trim_tails -h | --help

Options:
  -o --outfile=<STRING>  Name of output file
  -r --reverse           Trim beginning of sequence
  --penalty              Penalty against null  [default: 0]
  -t --target=<CHAR>     Target letter  [default: A]
  -h --help              Show this screen

"""

import numpy as np
from Bio import SeqIO

from util import new_record_seq_str


def _rev(s):
    return "".join(reversed(s))


def _likelihoods(seq, target, head_p, tail_p):
    head_lp = np.log(head_p)
    head_not_lp = np.log(1 - head_p)
    tail_lp = np.log(tail_p)
    tail_not_lp = np.log(1 - tail_p)

    bools = np.array(list(c == target for c in seq))
    head_counts = np.cumsum(bools)
    # need to insert an extra 0; since there is no head at idx=0
    head_counts = np.insert(head_counts, 0, 0)[:-1]
    tail_counts = np.cumsum(bools[::-1])[::-1]

    head_length = np.arange(len(seq))
    tail_length = head_length[::-1] + 1
    likelihoods = (head_counts * head_lp) + (head_length - head_counts) * head_not_lp \
        + (tail_counts * tail_lp) + (tail_length - tail_counts) * tail_not_lp

    n_target = bools.sum()
    null_log_likelihood = n_target * head_lp + (len(seq) - n_target) * head_not_lp

    return likelihoods, null_log_likelihood


def trim_tail(seq, target, head_p, tail_p, penalty=0, reverse=False):
    if reverse:
        seq = _rev(seq)
    likelihoods, null = _likelihoods(seq, target, head_p, tail_p)
    idx = np.argmax(likelihoods)
    max_log_likelihood = likelihoods[idx]
    result = seq
    tail = ""
    if (max_log_likelihood + penalty) > null:
        result, tail = seq[:idx], seq[idx:]
    if reverse:
        result = _rev(result)
        tail = _rev(tail)
    return result, tail


def trim_tails_file(infile, outfile, target, head_p, tail_p,
                    penalty=0, reverse=False, tailfile=None):
    records = list(SeqIO.parse(infile, 'fasta'))
    results = list(trim_tail(str(r.seq), target, head_p, tail_p,
                             penalty, reverse)
                   for r in records)

    keep, discard = zip(*results)
    keep_records = list(new_record_seq_str(r, s) for r, s in zip(records, keep))
    discard_records = list(new_record_seq_str(r, s) for r, s in zip(records, discard))

    # filter empty
    keep_records = list(r for r in keep_records if len(r))
    discard_records = list(r for r in discard_records if len(r))
    SeqIO.write(keep_records, outfile, 'fasta')
    if tailfile is not None:
        SeqIO.write(discard_records, tailfile, 'fasta')


if __name__ == "__main__":
    args = docopt(__doc__)
    filename = args["<infile>"]
    outfile = args["--outfile"]
    target = args["--target"]
    penalty = float(args["--penalty"])
    head_p = args["<head_p>"]
    tail_p = args["<tail_p>"]
    reverse = args["--reverse"]
    trim_tail_file(filename, outfile, target, head_p, tail_p,
                   penalty=penalty, reverse=reverse)
