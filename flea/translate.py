#!/usr/bin/env python
"""
Translate DNA reads from a fasta file.

"""

import sys

import click

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC, Gapped

from flea.util import insert_gaps


def _translate(record, gapped=False):
    result = record[:]
    if gapped:
        translated = record.seq.ungap('-').translate()
        result.seq = Seq(insert_gaps(str(record.seq), str(translated), '---', '-'),
                         alphabet=Gapped(IUPAC.IUPACProtein))
    else:
        result.seq = record.seq.translate()
    return result


def translate(infile, outfile, gapped=False):
    alphabet=IUPAC.ambiguous_dna
    if gapped:
        alphabet = Gapped(alphabet)
    records = SeqIO.parse(infile, "fasta", alphabet=alphabet)
    result = (_translate(r, gapped) for r in records)
    SeqIO.write(result, outfile, "fasta")


@click.command()
@click.option('-g', '--gapped', is_flag=True, help='allow gaps')
def main(gapped):
    translate(sys.stdin, sys.stdout, gapped)


if __name__ == "__main__":
    main()
