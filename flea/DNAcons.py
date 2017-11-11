#!/usr/bin/env python

"""Prints a consensus sequence.

Sequences must be aligned, and all the same length'

"""

import sys
import random
from collections import defaultdict

import click

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO.FastaIO import FastaWriter
from Bio import AlignIO
from Bio.Alphabet import IUPAC

from flea.util import grouper
from flea.util import id_to_copynumber


def _column_consensus(counter, seed=None):
    r = random
    if seed is not None:
        r = random.Random(int(seed))
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
    for seq, copynumber in zip(seqs, copies):
        for column, elt in enumerate(seq):
            counters[column][elt] += copynumber
    pairs = (_column_consensus(c, seed) for c in counters)
    cons, ambi = zip(*pairs)
    assert(len(cons) == len(ambi))
    for i, elt in enumerate(cons):
        assert(elt in ambi[i][1])
    cons = ''.join(cons)
    if len(cons) != aln_len:
        raise Exception('consensus has wrong length')
    return cons, ambi


@click.command()
@click.argument('filename')
@click.option('-o', '--outfile')
@click.option('--ambifile')
@click.option('--copynumbers', is_flag=True, help='seq id ends with copynumber')
@click.option('--codon', is_flag=True, help='do codon consensus.')
@click.option('--keep-gaps', is_flag=True, help='keep gaps in consensus')
@click.option('--name', help='record id for the fasta output')
@click.option('--seed', type=int, help='seed rng for breaking ties')
@click.option('-v', '--verbose', is_flag=True)
def consfile(filename, outfile=None, ambifile=None, copynumbers=False,
             codon=False, keep_gaps=False, name=None, verbose=False, seed=None):
    """Computes a consensus sequence and writes it to a file.

    Breaks ties independently for each position by choosing randomly
    among winners.

    If `ambifile` is not None, writes ambiguous calls to that file,
    one per line, in the following format:

    <0-index position> <frequency> <candidates>

    """
    ungap = not keep_gaps
    alignment = list(AlignIO.read(filename, "fasta"))
    seqs = list(r.seq for r in alignment)
    if copynumbers:
        copies = list(id_to_copynumber(r.id) for r in alignment)
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
    if name is None:
        name = "{}_cons".format(rec.name)
    if outfile is None:
        outfile = "{}.cons.fasta".format(filename)
    with open(outfile, 'w') as handle:
        writer = FastaWriter(handle, wrap=None)
        writer.write_header()
        newrecord = SeqRecord(Seq(str(_consensus), IUPAC.unambiguous_dna), id=name, description="")
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
        sys.stderr.write('Consensus written with name {}.\n'.format(name))


if __name__ == "__main__":
    consfile()
