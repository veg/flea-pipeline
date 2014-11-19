#!/usr/bin/env python
"""
Translate DNA reads from a fasta file.

Usage:
  translate.py <infile> <outfile>
  translate.py -h | --help

Options:
  -h --help     Show this screen.

"""

import sys

from docopt import docopt

from Bio import SeqIO
from Bio.Alphabet import IUPAC


def _translate(record):
    result = record[:]
    result.seq = record.seq.translate()
    return result


def translate(infile, outfile):
    records = SeqIO.parse(infile, "fasta", alphabet=IUPAC.ambiguous_dna)
    result = (_translate(r) for r in records)
    SeqIO.write(result, outfile, "fasta")


if __name__ == "__main__":
    args = docopt(__doc__)
    infile = args["<infile>"]
    outfile = args["<outfile>"]
    translate(infile, outfile)
