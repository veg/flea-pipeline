#!/usr/bin/env python

"""
Prints a gapless consensus sequence.

Sequences must be aligned, and all the same length'

Usage:
  DNAcons.py [options] <infile>
  DNAcons.py -h | --help

Options:
  -v --verbose            Print progress to STDERR
  --keep-gaps             Do not ungap the consensus
  --copynumbers=<STRING>  File containing "<id>\t<num>" lines.
  --id=<STRING>           Record id for the fasta output
  -o --outfile=<STRING>   Name of output file
  --ambifile=<STRING>
  -h --help               Show this screen

"""

import sys
import random
from collections import defaultdict

from docopt import docopt

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO.FastaIO import FastaWriter
from Bio import AlignIO
from Bio.Align import AlignInfo
from Bio.Alphabet import IUPAC


def _column_consensus(counter, seed=None):
    r = random
    if seed is not None:
        r = random.Random(seed)
    max_count = max(counter.values())
    cands = ''.join(sorted(char for char, count in counter.items()
                           if count == max_count))
    char = r.choice(cands)
    assert(char in cands)
    total = sum(counter.values())
    ambi = (max_count / total, cands)
    return char, ambi


def consensus(seqs, copies=None, seed=None):
    if copies is None:
        copies = [1] * len(seqs)
    counters = list(defaultdict(lambda: 0) for _ in range(len(seqs[0])))
    for seq, cn in zip(seqs, copies):
        for column, char in enumerate(seq):
            counters[column][char] += cn
    pairs = (_column_consensus(c, seed) for c in counters)
    cons, ambi = zip(*pairs)
    cons = ''.join(cons)
    assert(len(cons) == len(ambi))
    for i, char in enumerate(cons):
        assert(char in ambi[i][1])
    return cons, ambi


def consfile(filename, outfile=None, ambifile=None, copynumber_file=None,
             id_str=None, ungap=True, verbose=False, seed=None):
    """Computes a consensus sequence and writes it to a file.

    Breaks ties independently for each position by choosing randomly
    among winners.

    If `ambifile` is not None, writes ambiguous calls to that file,
    one per line, in the following format:

    <0-index position> <frequency> <candidates>

    """
    alignment = AlignIO.read(filename, "fasta")
    seqs = list(r.seq for r in alignment)
    if copynumber_file is None:
        copies = None
    else:
        with open(copynumber_file) as handle:
            lines = handle.read().strip().split('\n')
            pairs = list(e.split() for e in lines)
            cdict = dict((key, int(val)) for key, val in pairs)
            copies = list(cdict[r.id] for r in alignment)
    _consensus, ambiguous = consensus(seqs, copies, seed=seed)
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
    with open(outfile, 'w') as handle:
        writer = FastaWriter(handle, wrap=None)
        writer.write_header()
        newrecord = SeqRecord(Seq(str(_consensus), IUPAC.unambiguous_dna), id=id_str, description="")
        writer.write_record(newrecord)
    if ambifile is not None:
        with open(ambifile, 'w') as handle:
            for i, (freq, cands) in enumerate(ambiguous):
                if len(cands) > 1:
                    assert(_consensus[i] in cands)
                    handle.write('{} {} {}\n'.format(i, freq, cands))
    if verbose:
        sys.stderr.write('Consensus written with name {}.\n'.format(id_str))


if __name__ == "__main__":
    args = docopt(__doc__)
    filename = args["<infile>"]
    copynumber_file = args["--copynumbers"]
    verbose = args["--verbose"]
    id_str = args["--id"]
    outfile = args["--outfile"]
    ambifile = args['--ambifile']
    keep_gap = args["--keep-gaps"]
    ungap = not keep_gap
    consfile(filename, outfile, copynumber_file=copynumber_file,
             ambifile=ambifile, id_str=id_str, ungap=ungap, verbose=verbose)
