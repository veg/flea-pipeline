#!/usr/bin/env python

"""
Prints a gapless consensus sequence with any nucleotides lower that
50% identity recieving an ambiguous character X.

Sequences must be aligned, and all the same length'

Usage:
  DNAcons.py [options] <infile>
  DNAcons.py -h | --help

Options:
  -v --verbose           Print progress to STDERR
  --keep-gaps            Do not ungap the consensus
  --id=<STRING>          Record id for the fasta output
  -o --outfile=<STRING>  Name of output file
  --ambifile=<STRING>
  -h --help              Show this screen

"""

import sys
import random
from collections import Counter

from docopt import docopt

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO.FastaIO import FastaWriter
from Bio import AlignIO
from Bio.Align import AlignInfo
from Bio.Alphabet import IUPAC


def _column_consensus(counter, n_cols, seed=None):
    """

    >>> _column_consensus(Counter('aabc'), 4, seed=0)
    ('a', (0.5, ['a']))

    """
    r = random
    if seed is not None:
        r = random.Random(seed)
    mc = counter.most_common()
    count = mc[0][1]
    cands = sorted(c for c, n in mc if n == count)
    char = r.choice(cands)
    assert(char in cands)
    ambi = (count / n_cols, cands)
    return char, ambi


def consensus(seqs, seed=None):
    """

    >>> consensus(['aab', 'abb'], seed=0)[0]
    'abb'

    >>> consensus(['aab', 'abb'], seed=0)[1]
    ((1.0, ['a']), (0.5, ['a', 'b']), (1.0, ['b']))

    """
    pwm = list(Counter() for _ in range(len(seqs[0])))
    for seq in seqs:
        for i, char in enumerate(seq):
            pwm[i].update([char])
    n_cols = len(seqs)
    pairs = (_column_consensus(c, n_cols, seed) for c in pwm)
    cons, ambi = zip(*pairs)
    cons = ''.join(cons)
    assert(len(cons) == len(ambi))
    for i, char in enumerate(cons):
        assert(char in ambi[i][1])
    return ''.join(cons), ambi


def consfile(filename, outfile=None, ambifile=None, id_str=None,
             ungap=True, verbose=False):
    """Computes a consensus sequence and writes it to a file.

    Breaks ties independently for each position by choosing randomly
    among winners.

    If `ambifile` is not None, writes ambiguous calls to that file,
    one per line, in the following format:

    <0-index position> <frequency> <candidates>

    """
    alignment = AlignIO.read(filename, "fasta")
    seqs = list(r.seq for r in alignment)
    _consensus, ambiguous = consensus(seqs)
    if ungap:
        # cannot just call _consensus.ungap('-'), because that will
        # mess up ambiguity information
        indices = set(i for i, char in enumerate(_consensus) if char == '-')
        _consensus = ''.join(char for i, char in enumerate(_consensus)
                             if i not in indices)
        ambiguous = (pair for i, pair in enumerate(ambiguous)
                     if i not in indices)

    # use the first entry for the name, just so know what is what later
    rec = alignment[0]
    if id_str is None:
        id_str = "{}_cons".format(rec.name)
    if outfile is None:
        outfile = "{}.cons.fasta".format(filename)
    writer = FastaWriter(open(outfile, "w"), wrap=None)
    writer.write_header()
    newrecord = SeqRecord(Seq(str(_consensus), IUPAC.unambiguous_dna), id=id_str, description="")
    writer.write_record(newrecord)
    if ambifile is not None:
        with open(ambifile, 'w') as handle:
            for i, (freq, cands) in enumerate(ambiguous):
                if len(cands) > 1:
                    assert(_consensus[i] in cands)
                    handle.write('{} {} {}\n'.format(i, freq, ''.join(cands)))
    if verbose:
        sys.stderr.write('Consensus written with name {}.\n'.format(id_str))


if __name__ == "__main__":
    args = docopt(__doc__)
    filename = args["<infile>"]
    verbose = args["--verbose"]
    id_str = args["--id"]
    outfile = args["--outfile"]
    ambifile = args['--ambifile']
    keep_gap = args["--keep-gaps"]
    ungap = not keep_gap
    consfile(filename, outfile, ambifile, id_str, ungap=ungap, verbose=verbose)

