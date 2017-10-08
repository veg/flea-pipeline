#!/usr/bin/env python

"""
Trim poly-A and poly-T heads and tails. Uses a simple HMM.

Supports both FASTA and FASTQ files. Default is FASTA; use --fastq
flag to switch.

Usage:
  trim_tails.py [options] <infile> <outfile>
  trim_tails.py -h | --help

Options:
  --fastq  Treat input and output files as FASTQ.
  --train  Run Baum-Welch training.
  --jobs <INT>  Run in parallel.  [default: 1]
  --max-iters <INT>  Max EM iterations.  [default: 16]
  -h --help

"""
import os
import multiprocessing

from docopt import docopt
import numpy as np
from Bio import SeqIO
from pomegranate import DiscreteDistribution, State, HiddenMarkovModel

from flea.util import new_record_seq_str

# states
P = 0.985
Q = 0.005
head_tail_dist_A = DiscreteDistribution({'A': P, 'C': Q, 'G': Q, 'T': Q})
head_tail_dist_T = DiscreteDistribution({'A': Q, 'C': Q, 'G': Q, 'T': P})

# freeze the body distribution
body_dist = DiscreteDistribution({'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25}, frozen=True)

# heads and tails get same Python object for their distribution, which
# ties them together.
# TODO: how to also tie the P and Q values inside each distribution?
HEAD_A = State(head_tail_dist_A, name='head_a')
HEAD_T = State(head_tail_dist_T, name='head_t')
BODY = State(body_dist, name='body')
TAIL_A = State(head_tail_dist_A, name='tail_a')
TAIL_T = State(head_tail_dist_T, name='tail_t')

# model
hmm = HiddenMarkovModel()
hmm.add_states(HEAD_A, HEAD_T, BODY, TAIL_A, TAIL_T)

# start probs
hmm.add_transition(hmm.start, HEAD_A, 0.1, group='a')
hmm.add_transition(hmm.start, HEAD_T, 0.1, group='a')
hmm.add_transition(hmm.start, BODY, 0.8)

# transitions
# tie head and tail transitions together
hmm.add_transition(HEAD_A, HEAD_A, 0.99)
hmm.add_transition(HEAD_T, HEAD_T, 0.99)
hmm.add_transition(HEAD_A, BODY, 0.01, group='b')
hmm.add_transition(HEAD_T, BODY, 0.01, group='b')
hmm.add_transition(BODY, BODY, 0.99)
hmm.add_transition(BODY, TAIL_A, 0.005, group='c')
hmm.add_transition(BODY, TAIL_T, 0.005, group='c')
hmm.add_transition(TAIL_A, TAIL_A, 1.0)
hmm.add_transition(TAIL_T, TAIL_T, 1.0)

hmm.bake()


def trim(seq):
    """Split into head, body, and tail.

    >>> trim("AAAAAAGTACTTTTT")
    ('AAAAAA', 'GTAC', 'TTTTT')

    >>> trim("TTCTTCCCA")
    ('', 'TTCTTCCCA', '')

    """
    states = hmm.predict(np.array(list(seq)))
    head, body, tail = [], [], []
    for s_id, c in zip(states, seq):
        s = hmm.states[s_id]
        if s in [HEAD_A, HEAD_T]:
            head.append(c)
        elif s == BODY:
            body.append(c)
        elif s in [TAIL_A, TAIL_T]:
            tail.append(c)
        else:
            raise Exception('unknown state: {}'.format(s))
    return ''.join(head), ''.join(body), ''.join(tail)


def trim_record(r):
    fastq = 'phred_quality' in r.letter_annotations
    head, body, tail = trim(str(r.seq))

    head_record = new_record_seq_str(r, head)
    body_record = new_record_seq_str(r, body)
    tail_record = new_record_seq_str(r, tail)

    if fastq:
        quality = r.letter_annotations['phred_quality']

        head_record.letter_annotations['phred_quality'] = quality[:len(head)]
        body_record.letter_annotations['phred_quality'] = quality[len(head):(len(head) + len(body))]
        tail_record.letter_annotations['phred_quality'] = quality[(len(head) + len(body)):]
    return head_record, body_record, tail_record


def _notempty(rs):
    return (r for r in rs if len(r.seq))


def trim_file(infile, outfile, fastq=False, train=False, max_iters=16, n_jobs=1):
    filetype = 'fastq' if fastq else 'fasta'
    records = list(SeqIO.parse(infile, filetype))

    # train model
    if train:
        seqs = list(str(rec.seq).upper() for rec in records)
        hmm.fit(seqs, max_iterations=16, verbose=False, n_jobs=n_jobs)

    if n_jobs > 1:
        pool = multiprocessing.Pool(n_jobs)
        results = pool.map(trim_record, records)
    else:
        results = list(trim_record(rec) for rec in records)

    heads, bodies, tails = zip(*results)

    basename, ext = os.path.splitext(outfile)
    head_file = "{}.heads{}".format(basename, ext)
    tail_file = "{}.tails{}".format(basename, ext)
    SeqIO.write(_notempty(heads), head_file, filetype)
    SeqIO.write(_notempty(bodies), outfile, filetype)
    SeqIO.write(_notempty(tails), tail_file, filetype)

    model_json = hmm.to_json()
    with open("{}.hmm.json".format(basename), 'w') as handle:
        handle.write(model_json)


if __name__ == "__main__":
    args = docopt(__doc__)
    filename = args["<infile>"]
    outfile = args["<outfile>"]
    fastq = args['--fastq']
    train = args['--train']
    max_iters = int(args['--max-iters'])
    n_jobs = int(args['--jobs'])
    trim_file(filename, outfile, fastq=fastq, train=train,
              max_iters=max_iters, n_jobs=n_jobs)
