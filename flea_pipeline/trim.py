#!/usr/bin/env python

"""
Trim poly-A and poly-T heads and tails. Uses a simple HMM.

Usage:
  trim_tails.py <infile> <outfile>
  trim_tails.py -h | --help

Options:
  -h --help

"""
import warnings

import numpy as np
from Bio import SeqIO
from docopt import docopt

from flea_pipeline.util import new_record_seq_str

# FIXME: do not hardcode these values

# STATES
A_HEAD = 0
T_HEAD = 1
BODY = 2
A_TAIL = 3
T_TAIL = 4

HEAD_SELF = 0.99
BODY_SELF = 0.99
TO_TAIL = (1 - BODY_SELF) / 2
TAIL_SELF = 1
TRANSMAT = np.array([[HEAD_SELF, 0, 1 - HEAD_SELF, 0, 0],
                     [0, HEAD_SELF, 1 - HEAD_SELF, 0, 0],
                     [0, 0, BODY_SELF, TO_TAIL, TO_TAIL],
                     [0, 0, 0, TAIL_SELF, 0],
                     [0, 0, 0, 0, TAIL_SELF]])

P_NO_HEAD = 0.8
P_HEAD = (1 - P_NO_HEAD) / 2
STARTPROB = np.array([P_HEAD, P_HEAD, P_NO_HEAD, 0, 0])

# SYMBOLS
# 0 : A
# 1 : C
# 2 : G
# 3 : T
P = 0.985
Q = 0.005
E = 0.25
EMISSIONPROB = np.array([[P, Q, Q, Q],
                         [Q, Q, Q, P],
                         [E, E, E, E],
                         [P, Q, Q, Q],
                         [Q, Q, Q, P]])



def decode(obs, transmat, emissionprob, startprob):
    """Viterbi algorithm.

    transmat: (n_states x n_states) [i, j]: transition from i to j
    emissionprob: (n_states x n_symbols)
    startprob: (n_states,)

    """
    startprob = startprob.ravel()
    n_states = transmat.shape[0]
    n_symbols = emissionprob.shape[1]
    if transmat.shape != (n_states, n_states):
        raise Exception('Transmission matrix shape {} is not square.'.format(transmat.shape))
    if len(startprob) != n_states:
        raise Exception('Wrong number of starting probabilities.'
                        ' Expected {} and got {}.'.format(n_states, len(startprob)))
    if emissionprob.shape != (n_states, n_symbols):
        raise Exception('Emission matrix has wrong number of states.'
                        ' Expected {} and got {}.'.format(n_states, emissionprob.shape[0]))

    if not np.all(np.sum(transmat, axis=1) == 1):
        raise Exception('Transmission probabilities do not sum to 1.')
    if not np.sum(startprob) == 1:
        raise Exception('Starting probabilities do not sum to 1.')
    if not np.all(np.sum(emissionprob, axis=1) == 1):
        raise Exception('Emission probabilities do not sum to 1.')

    if np.max(obs) > n_symbols:
        raise Exception('Observation contains an invalid state: {}'.format(np.max(obs)))

    with warnings.catch_warnings():
        # already checked probabilities, so should be safe to ignoring warnings
        warnings.simplefilter("ignore")
        logtrans = np.log(transmat)
        logstart = np.log(startprob)
        logemission = np.log(emissionprob)

    # heavily adapted from:
    # http://phvu.net/2013/12/06/sweet-implementation-of-viterbi-in-python/

    # dynamic programming matrix
    trellis = np.zeros((n_states, len(obs)))
    # back pointers
    backpt = np.ones((n_states, len(obs)), np.int) * -1

    trellis[:, 0] = logstart + logemission[:, obs[0]]
    for t in range(1, len(obs)):
        logprobs = trellis[:, t - 1].reshape(-1, 1) + logtrans + logemission[:, obs[t]].reshape(1, -1)
        trellis[:, t] = logprobs.max(axis=0)
        backpt[:, t] = logprobs.argmax(axis=0)
    result = [trellis[:, -1].argmax()]
    for i in range(len(obs)-1, 0, -1):
        result.append(backpt[result[-1], i])
    result = result[::-1]
    return result


def trim(seq):
    """Split into head, body, and tail.

    >>> trim("AAAAAGTACTTTTT")
    ('AAAAA', 'GTAC', 'TTTTT')

    >>> trim("TTCTTCCCA")
    ('', 'TTCTTCCCA', '')

    """
    nuc_map = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    obs = np.array(list(nuc_map[c.upper()] for c in seq))
    states = decode(obs, TRANSMAT, EMISSIONPROB, STARTPROB)
    head, body, tail = [], [], []
    for s, c in zip(states, seq):
        if s in [A_HEAD, T_HEAD]:
            head.append(c)
        elif s == BODY:
            body.append(c)
        elif s in [A_TAIL, T_TAIL]:
            tail.append(c)
        else:
            raise Exception('unknown state: {}'.format(s))
    return ''.join(head), ''.join(body), ''.join(tail)


def _notempty(g):
    return (x for x in g if len(x))


def trim_file(infile, outfile):
    records = list(SeqIO.parse(infile, 'fasta'))
    results = list(trim(str(r.seq)) for r in records)
    head, body, tail = zip(*results)

    head_records = _notempty(new_record_seq_str(r, s) for r, s in zip(records, head))
    body_records = _notempty(new_record_seq_str(r, s) for r, s in zip(records, body))
    tail_records = _notempty(new_record_seq_str(r, s) for r, s in zip(records, tail))

    head_file = "{}.heads".format(outfile)
    tail_file = "{}.tails".format(outfile)
    SeqIO.write(body_records, outfile, 'fasta')
    SeqIO.write(head_records, head_file, 'fasta')
    SeqIO.write(tail_records, tail_file, 'fasta')


if __name__ == "__main__":
    args = docopt(__doc__)
    filename = args["<infile>"]
    outfile = args["<outfile>"]
    trim_file(filename, outfile)
