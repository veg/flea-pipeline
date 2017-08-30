#!/usr/bin/env python

"""
Prints a consensus sequence.

Sequences must be aligned, and all the same length'

Usage:
  DNAcons.py [options] <infile>
  DNAcons.py -h | --help

Options:
  -v --verbose            Print progress to STDERR
  --keep-gaps             Do not ungap the consensus
  --copynumbers           Seq id ends with copynumber
  --id=<STRING>           Record id for the fasta output
  -o --outfile=<STRING>   Name of output file
  --ambifile=<STRING>
  --codon                 Do codon alignment.
  -h --help               Show this screen

"""

import sys
import random
from collections import defaultdict

from docopt import docopt

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO.FastaIO import FastaWriter
from Bio import AlignIO
from Bio.Alphabet import IUPAC

from flea.util import grouper
from flea.util import id_to_con


def _column_consensus(counter, seed=None, codon=False):
    r = random
    if seed is not None:
        r = random.Random(seed)
    max_count = max(counter.values())
    cands = list(sorted(elt for elt, count in counter.items()
                        if count == max_count))
    elt = r.choice(cands)
    assert(elt in cands)
    total = sum(counter.values())
    ambi = (max_count / total, sorted(cands))
    return elt, ambi


def consensus(seqs, copies=None, codon=False, seed=None):
    aln_len = len(seqs[0])
    if not all(len(s) == aln_len for s in seqs):
        raise Exception('sequences are not aligned')
    if codon:
        if not aln_len % 3 == 0:
            raise Exception('alignment length is not multiple of 3')
        seqs = list(list("".join(g) for g in grouper(s, 3)) for s in seqs)
    if copies is None:
        copies = [1] * len(seqs)
    counters = list(defaultdict(lambda: 0) for _ in range(len(seqs[0])))
    for seq, cn in zip(seqs, copies):
        for column, elt in enumerate(seq):
            counters[column][elt] += cn
    pairs = (_column_consensus(c, seed, codon) for c in counters)
    cons, ambi = zip(*pairs)
    assert(len(cons) == len(ambi))
    for i, elt in enumerate(cons):
        assert(elt in ambi[i][1])
    cons = ''.join(cons)
    return cons, ambi


def consfile(filename, outfile=None, ambifile=None, do_copynumber=False,
             id_str=None, codon=False, ungap=True, verbose=False, seed=None):
    """Computes a consensus sequence and writes it to a file.

    Breaks ties independently for each position by choosing randomly
    among winners.

    If `ambifile` is not None, writes ambiguous calls to that file,
    one per line, in the following format:

    <0-index position> <frequency> <candidates>

    """
    alignment = list(AlignIO.read(filename, "fasta"))
    seqs = list(r.seq for r in alignment)
    if do_copynumber:
        copies = list(id_to_cn(r.id) for r in alignment)
    else:
        copies = None
    _consensus, ambiguous = consensus(seqs, copies, codon=codon, seed=seed)
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
                    if codon:
                        chosen = _consensus[i * 3: i * 3 + 3]
                    else:
                        chosen = _consensus[i]
                    assert(chosen in cands)
                    handle.write('{} {} {}\n'.format(i, freq, cands))
    if verbose:
        sys.stderr.write('Consensus written with name {}.\n'.format(id_str))


if __name__ == "__main__":
    args = docopt(__doc__)
    filename = args["<infile>"]
    do_copynumber = args["--copynumbers"]
    verbose = args["--verbose"]
    id_str = args["--id"]
    outfile = args["--outfile"]
    ambifile = args['--ambifile']
    keep_gap = args["--keep-gaps"]
    codon = args["--codon"]
    ungap = not keep_gap
    consfile(filename, outfile, do_copynumber=do_copynumber,
             ambifile=ambifile, id_str=id_str, codon=codon, ungap=ungap,
             verbose=verbose)
