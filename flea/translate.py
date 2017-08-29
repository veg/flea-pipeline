#!/usr/bin/env python
"""
Translate DNA reads from a fasta file.

Usage:
  translate.py [options] < infile  > outfile
  translate.py -h | --help

Options:
  -g --gapped   Allow gaps.  [default: False]
  -h --help     Show this screen.

"""

import sys
import argparse

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


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="My parser")

    parser.add_argument('--gapped', dest='gapped', action='store_true')
    parser.add_argument('--no-gapped', dest='gapped', action='store_false')
    parser.set_defaults(feature=False)

    parsed = parser.parse_args()
    
    translate(sys.stdin, sys.stdout, parsed.gapped)
