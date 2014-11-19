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


if __name__ == "__main__":
    args = docopt(__doc__)
    infile = args["<infile>"]
    outfile = args["<outfile>"]
    records = SeqIO.parse(infile, "fasta", alphabet=IUPAC.unambiguous_dna)
    result = (record.translate() for r in records)
    SeqIO.write(records, outfile, "fasta")
